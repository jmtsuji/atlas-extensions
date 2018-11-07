#!/usr/bin/env bash
set -euo pipefail

# If no input is provided, provide help and exit
if [ $# == 0 ]; then
	# Assign script name
	script_name=${0##*/}
	script_name=${0%.*}

	# Help statement
	printf "${0##*/}: Wrapper to calculate the read recruitment of genome bins to metagenomes.\n"
	printf "Contact Jackson M. Tsuji (jackson.tsuji@uwaterloo.ca) for bug reports or feature requests.\n\n"
	printf "Installation: you can install all dependencies via conda and then add the scripts from the Github repo with:\n"
	printf "              conda create -y -n bin_mapping_stats -c conda-forge -c bioconda -c r samtools bbmap r r-plyr r-dplyr r-getopt\n"
	printf "              source activate bin_mapping_stats\n"
	printf "              git clone https://github.com/jmtsuji/atlas-extensions.git\n"
	printf "              cd atlas-extensions\n"
	printf "              git checkout perfectmode\n" # TODO - eliminate this as the script is merged into the master
	printf "              env_path=\$(which samtools)\n"
	echo "              env_path=\${env_path%/*}"
	printf "              cp calculate_bin_abundance_in_metagenome.sh calculate_coverage_stats.R aggregate_mapping_stats.R \${env_path}\n"
	printf "Usage: ${0##*/} refined_bin_dir raw_read_dir output_dir threads memory 2>&1 | tee ${script_name}.log\n\n"
	printf "Usage details:\n"
	printf "   1. refined_bin_dir: Directory containing genomes/bins. MUST have extension *.fa!! Symlinks are okay.\n"
	printf "   2. raw_read_dir: Directory containing QC'ed, unassembled metagenomic reads with naming structure matching that of the ATLAS pipeline.\n"
	printf "                    *QC_R1.fastq.gz, *QC_R2.fastq.gz, and *QC_se.fastq.gz all needed. Symlinks are okay.\n"
	printf "   3. output_dir: Directory where you want the results to be output. Anything already there might be overwritten.\n"
	printf "   4. threads: number of threads to run. Caps out for performance increase at around 12-16.\n"
	printf "   5. memory: in gigabytes.\n\n"

	# Exit
	exit 1
fi

# Set user variables
refined_bin_dir=$1 # bins MUST have extension *.fa
raw_read_dir=$2 # reads MUST be output of ATLAS QC and match naming structure. R1, R2, se reads all need to be there, for all samples.
output_dir=$3
threads=$4
memory=$5 # in Gigabytes

# Startup reporting
(>&2 echo "[ $(date -u) ]: Running ${0##*/}")
(>&2 echo "[ $(date -u) ]: Command: ${0##*/} ${@}")
(>&2 echo "[ $(date -u) ]: refined_bin_dir: ${refined_bin_dir}")
(>&2 echo "[ $(date -u) ]: raw_read_dir: ${raw_read_dir}")
(>&2 echo "[ $(date -u) ]: output_dir: ${output_dir}")
(>&2 echo "[ $(date -u) ]: bin_mapping_summary_filename: bin_mapping_stats.tsv (in output_dir)")
(>&2 echo "[ $(date -u) ]: threads: ${threads}")
(>&2 echo "[ $(date -u) ]: memory: ${memory} GB")

###### Script setup ######
# Create folder structure in output dir
mkdir -p ${output_dir}/mapping/logs \
	${output_dir}/mapping/lists/R1 \
	${output_dir}/mapping/lists/R2 \
	${output_dir}/mapping/bam \
	${output_dir}/coverage/by_nucleotide \
	${output_dir}/coverage/by_contig \
	${output_dir}/coverage/logs \
	${output_dir}/detailed_stats

# Assign names to stats summary tables
genome_stats_filename="${output_dir}/detailed_stats/genome_stats.tsv"
read_count_summary_filename="${output_dir}/detailed_stats/metagenome_read_counts.tsv"
bin_mapping_summary_filename="${output_dir}/detailed_stats/bin_mapping_stats.tsv"
coverage_summary_filename="${output_dir}/detailed_stats/coverage_summary.tsv"
final_stats_summary="${output_dir}/genome_bin_mapping_stats.tsv"
mapped_read_list_summary_filename="${output_dir}/detailed_stats/mapped_read_list_metadata.tsv"

# Find bin and raw read files
bin_paths=($(find -L ${refined_bin_dir} -iname "*.fa")) # Careful - make sure the files match *.fa!
raw_read_paths=($(find -L ${raw_read_dir} -iname "*R1.fastq.gz")) # Grab the R1 only
iterations=$((${#bin_paths[@]}*${#raw_read_paths[@]}))

(>&2 echo "[ $(date -u) ]: found ${#bin_paths[@]} genome bins with extension *.fa in the refined_bin_dir")
(>&2 echo "[ $(date -u) ]: found ${#raw_read_paths[@]} set of unassembled metagenomic reads in the raw_read_dir")


###### Part 1: calculate bin stats ######
(>&2 echo "[ $(date -u) ]: collecting genome bin statistics using statswrapper.sh (-> ${genome_stats_filename##*/})")
(>&2 printf "[ $(date -u) ]: ")
statswrapper.sh ${refined_bin_dir}/*.fa format=5 > ${genome_stats_filename}


###### Part 2: count metagenome reads ######
# Initialize output metagenome read count file
(>&2 echo "[ $(date -u) ]: Initiatilizing '${read_count_summary_filename##*/}'")
printf "metagenome\ttotal_reads\tR1_reads\tR2_reads\tse_reads\n" > ${read_count_summary_filename}

# Get a count of the number of total reads in each metagenome
(>&2 echo "[ $(date -u) ]: Counting total reads in the ${#raw_read_paths[@]} provided metagenomes")
for raw_read_path in ${raw_read_paths[@]}; do

	# Get base names of reads
	raw_read_name_base=${raw_read_path%*_QC_R1.fastq.gz}
	raw_read_name_base=${raw_read_name_base##*/}

	# Assign generic variables for read mapping
	R1=${raw_read_dir}/${raw_read_name_base}_QC_R1.fastq.gz
	R2=${raw_read_dir}/${raw_read_name_base}_QC_R2.fastq.gz
	se=${raw_read_dir}/${raw_read_name_base}_QC_se.fastq.gz

	(>&2 printf "[ $(date -u) ]: Counting ${raw_read_name_base}*.fastq.gz: ")

	# Count total reads (rough)
	R1_count=$(($(zcat ${R1} | wc -l)/4))
	R2_count=$(($(zcat ${R2} | wc -l)/4))
	se_count=$(($(zcat ${se} | wc -l)/4))
	total_count=$((${R1_count}+${R2_count}+${se_count}))
	(>&2 echo "${R1_count} + ${R2_count} + ${se_count} = ${total_count} reads")

	# Store in data frame
	printf "${raw_read_name_base}\t${total_count}\t${R1_count}\t${R2_count}\t${se_count}\n" >> ${read_count_summary_filename}

done
(>&2 echo "[ $(date -u) ]: metagenome read counting complete")


###### Part 3: map reads ######
# Initialize output read stats file
(>&2 echo "[ $(date -u) ]: Initiatilizing '${bin_mapping_summary_filename##*/}'")
printf "genome\tmetagenome\tmapped_reads\n" > ${bin_mapping_summary_filename}

# Initialize coverage summary file
echo "[ $(date -u) ]: Initializing '${coverage_summary_filename}'"
printf "genome\tmetagenome\tcoverage_mean\tcoverage_sd\tpercent_contigs_with_zero_coverage_event\tcoverage_mean_filtered\tcoverage_sd_filtered\tpercent_contigs_with_zero_coverage_event_filtered\n" > ${coverage_summary_filename}

# Inititalize Emilie table
echo "[ $(date -u) ]: Initializing '${mapped_read_list_summary_filename}'"
printf "genome\tmetagenome\tmapped_read_list_filepath_R1\tmapped_read_list_filepath_R2\n" > ${mapped_read_list_summary_filename}

# Start counting the number of iterations processed
iteration=1

# Run all against all in a nested 'for' loop
(>&2 echo "[ $(date -u) ]: Mapping metagenome reads to genomes for ${#bin_paths[@]} * ${#raw_read_paths[@]} = ${iterations} combinations")
for bin_path in ${bin_paths[@]}; do

	# Get base name of bin
	bin_name_base=${bin_path%.*}
	bin_name_base=${bin_name_base##*/}

	for raw_read_path in ${raw_read_paths[@]}; do

		# Get base names of reads
		raw_read_name_base=${raw_read_path%*_QC_R1.fastq.gz}
		raw_read_name_base=${raw_read_name_base##*/}

		# Assign generic variables for read mapping
		R1=${raw_read_dir}/${raw_read_name_base}_QC_R1.fastq.gz
		R2=${raw_read_dir}/${raw_read_name_base}_QC_R2.fastq.gz
		se=${raw_read_dir}/${raw_read_name_base}_QC_se.fastq.gz
		samtools_depth_filename="${output_dir}/coverage/by_nucleotide/${raw_read_name_base}_to_${bin_name_base}.tsv"
		mapping_logfile=${output_dir}/mapping/logs/${raw_read_name_base}_to_${bin_name_base}_contig_coverage_stats.log
		bam_filename_base="${raw_read_name_base}_to_${bin_name_base}"
		bam_filename="${output_dir}/mapping/bam/${bam_filename_base}.bam"

		# Read map AND pipe directly to stats (to avoid excessive input/output, which is rough on hard drives)
		(>&2 printf "[ $(date -u) ]: ${iteration}: mapping '${raw_read_name_base}*fastq.gz' to '${bin_name_base}': ")
		bbwrap.sh nodisk=t ref=${bin_path} in1=${R1},${se} in2=${R2},null perfectmode=t trimreaddescriptions=t \
			outm=stdout threads=${threads} pairlen=1000 pairedonly=t mdtag=t xstag=fs nmtag=t sam=1.3 \
			local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${memory}G 2> ${mapping_logfile} | \
			samtools view -@ ${threads} -O bam | samtools sort -@ ${threads} 2>/dev/null > ${bam_filename}

		# Calculate the number of mapped reads
		mapped_reads=$(samtools view -c -F 4 ${bam_filename})
		(>&2 echo "${mapped_reads} mapped reads")

		# Add to TSV file
		printf "${bin_name_base}\t${raw_read_name_base}\t${mapped_reads}\n" >> ${bin_mapping_summary_filename}
		
		# Generate coverage profile
		samtools depth -aa ${bam_filename} > ${samtools_depth_filename}

		# For Emilie: export mapped reads and make data table
		# Keep R1 and R2 separate to avoid pseudoreplication
		mapped_read_list_dir_R1="${output_dir}/mapping/lists/R1/${raw_read_name_base}" # TODO - set these variables earlier in the code for clarity
		mapped_read_list_dir_R2="${output_dir}/mapping/lists/R2/${raw_read_name_base}" # TODO - set these variables earlier in the code for clarity
		mapped_read_list_filepath_R1="${mapped_read_list_dir_R1}/${bin_name_base}.R1.list"
		mapped_read_list_filepath_R2="${mapped_read_list_dir_R2}/${bin_name_base}.R2.list"
		(>&2 printf "[ $(date -u) ]: ${iteration}: Generating summary of mapped reads\n")
		mkdir -p ${mapped_read_list_dir} # TODO - clean this up
		samtools view -F 4 -f 64 ${bam_filename} | cut -d $'\t' -f 1 > ${mapped_read_list_filepath_R1}
		samtools view -F 4 -f 128 ${bam_filename} | cut -d $'\t' -f 1 > ${mapped_read_list_filepath_R2}
		# Add metadata to a TSV table
		printf "${bin_name_base}\t${raw_read_name_base}\t${mapped_read_list_filepath_R1}\t${mapped_read_list_filepath_R2}\n" >> ${mapped_read_list_summary_filename}


		# Only determine coverage stats if mapped reads > 0
		if [ ${mapped_reads} -gt 0 ]; then

			# Set some variables for coverage stats
			coverage_by_contig_filename="${output_dir}/coverage/by_contig/${raw_read_name_base}_to_${bin_name_base}.tsv"
			coverage_logfile="${output_dir}/coverage/logs/${raw_read_name_base}_to_${bin_name_base}_summary.log"
			zero_cov_threshold=100 # TODO - HARD-CODED for now!

			# Roughly estimate read length of metagenome
			# TODO - improve this.
			# Temporarily disable script exiting with non-normal exit statuses. 'Head' may be causing problems - see https://stackoverflow.com/a/19120674 (accessed 180927)
			set +e
			read_length=$(($(zcat ${R1} | head -n 2 | tail -n 1 | wc -m)-1))
			set -e

			# Summarize coverage stats
			(>&2 echo "[ $(date -u) ]: ${iteration}: summarizing coverage stats (assumed read length of ${read_length})")
			calculate_coverage_stats.R --samtools_coverage_table ${samtools_depth_filename} --bin_ID ${bin_name_base} \
				--metagenome_ID ${raw_read_name_base} --output_contig_stats_table ${coverage_by_contig_filename} \
				--read_length ${read_length} --zero_coverage_threshold ${zero_cov_threshold} 2> ${coverage_logfile} | \
				tail -n 1 >> ${coverage_summary_filename}

		else

			(>&2 echo "[ $(date -u) ]: ${iteration}: Skipping coverage stats summary due to 0 mapped reads")
			# Make empty output files
			touch "${output_dir}/coverage/by_contig/${raw_read_name_base}_to_${bin_name_base}.tsv" \
				"${output_dir}/coverage/logs/${raw_read_name_base}_to_${bin_name_base}_summary.log"

		fi

		# Increase the iterator for tracking the number of processed genome-metagenome pairs
		iteration_tmp=$((${iteration}+1))
		iteration=${iteration_tmp}

	done
done
(>&2 echo "[ $(date -u) ]: read mapping complete")


###### Part 4: aggregate stats ######
(>&2 echo "[ $(date -u) ]: aggregating stats (-> ${final_stats_summary##*/})")
aggregate_mapping_stats.R -g ${genome_stats_filename} \
	-r ${read_count_summary_filename} \
	-m ${bin_mapping_summary_filename} \
	-c ${coverage_summary_filename} \
	-o ${final_stats_summary} 2>/dev/null


(>&2 echo "[ $(date -u) ]: ${0##*/}: Finished.")


