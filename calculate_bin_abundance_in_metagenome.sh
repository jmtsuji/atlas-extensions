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
	printf "              conda create -y -n calc_bin_abundance -c conda-forge -c bioconda samtools bbmap\n"
	printf "              source activate calc_bin_abundance\n\n"
	printf "Usage: ${0##*/} refined_bin_dir raw_read_dir output_dir threads memory 2>&1 | tee ${script_name}.log\n\n"
	printf "Usage details:\n"
	printf "   1. refined_bin_dir: Directory containing genomes/bins. MUST have extension *.fa!! Symlinks are okay.\n"
	printf "   2. raw_read_dir: Directory containing QC'ed, unassembled metagenomic reads with naming structure matching that of the ATLAS pipeline.\n"
	printf "                    *QC_R1.fastq.gz, *QC_R2.fastq.gz, and *QC_se.fastq.gz all needed. Symlinks are okay.\n"
	printf "   3. output_dir: Directory where you want the results to be output. Anything already there might be overwritten.\n"
	printf "   4. threads: number of threads to run.\n"
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
(>&2 echo "[ $(date -u) ]: refined_bin_dir: ${refined_bin_dir}")
(>&2 echo "[ $(date -u) ]: raw_read_dir: ${raw_read_dir}")
(>&2 echo "[ $(date -u) ]: output_dir: ${output_dir}")
(>&2 echo "[ $(date -u) ]: bin_mapping_summary_filename: bin_mapping_stats.tsv (in output_dir)")
(>&2 echo "[ $(date -u) ]: threads: ${threads}")
(>&2 echo "[ $(date -u) ]: memory: ${memory} GB")

# Create folder structure in output dir
mkdir -p ${output_dir}/mapping ${output_dir}/logs ${output_dir}/coverage ${output_dir}/detailed_stats
read_count_summary_filename="${output_dir}/detailed_stats/metagenome_read_counts.tsv"
bin_mapping_summary_filename="${output_dir}/detailed_stats/bin_mapping_stats.tsv"

# Find bin and raw read files
bin_paths=($(find -L ${refined_bin_dir} -iname "*.fa")) # Careful - make sure the files match *.fa!
raw_read_paths=($(find -L ${raw_read_dir} -iname "*R1.fastq.gz")) # Grab the R1 only
iterations=$((${#bin_paths[@]}*${#raw_read_paths[@]}))

(>&2 echo "[ $(date -u) ]: found ${#bin_paths[@]} genome bins with extension *.fa in the refined_bin_dir")
(>&2 echo "[ $(date -u) ]: found ${#raw_read_paths[@]} set of unassembled metagenomic reads in the raw_read_dir")

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

# Initialize output read stats file
(>&2 echo "[ $(date -u) ]: Initiatilizing '${bin_mapping_summary_filename##*/}'")
printf "genome\tmetagenome\tmapped_reads\tgenome_length_nt\n" > ${bin_mapping_summary_filename}

# Start counting the number of iterations processed
iteration=1

# Run all against all in a nested 'for' loop
(>&2 echo "[ $(date -u) ]: Mapping metagenome reads to genomes for ${#bin_paths[@]} * ${#raw_read_paths[@]} = ${iterations} combinations")
for bin_path in ${bin_paths[@]}; do

	# Get base name of bin
	bin_name_base=${bin_path%.*}
	bin_name_base=${bin_name_base##*/}

	# Get total nucleotide length in bin for later
	genome_char_count=$(grep -v "^>" ${bin_path} | wc)
	genome_new_lines=$(echo ${genome_char_count} | cut -d ' ' -f 1)
	genome_total_chars=$(echo ${genome_char_count} | cut -d ' ' -f 3)
	genome_length_nt=$((${genome_total_chars}-${genome_new_lines}))

	(>&2 echo "[ $(date -u) ]: Starting on bin '${bin_path##*/}' with length ${genome_length_nt} nt")

	for raw_read_path in ${raw_read_paths[@]}; do

		# Get base names of reads
		raw_read_name_base=${raw_read_path%*_QC_R1.fastq.gz}
		raw_read_name_base=${raw_read_name_base##*/}

		# Assign generic variables for read mapping
		R1=${raw_read_dir}/${raw_read_name_base}_QC_R1.fastq.gz
		R2=${raw_read_dir}/${raw_read_name_base}_QC_R2.fastq.gz
		se=${raw_read_dir}/${raw_read_name_base}_QC_se.fastq.gz
		logfile=${output_dir}/logs/${bin_name_base}_to_${raw_read_name_base}_contig_coverage_stats.log

		# Read map AND pipe directly to stats (to avoid excessive input/output, which is rough on hard drives)
		(>&2 printf "[ $(date -u) ]: ${iteration}: mapping '${raw_read_name_base}*fastq.gz' to '${bin_name_base}': ")
		bbwrap.sh nodisk=t ref=${bin_path} in1=${R1},${se} in2=${R2},null perfectmode=t trimreaddescriptions=t \
			out=stdout threads=${threads} pairlen=1000 pairedonly=t mdtag=t xstag=fs nmtag=t sam=1.3 \
			local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${memory}G 2> ${logfile} | \
			tee >(samtools view -@ $((${threads}/2)) -c -F 4 > ${output_dir}/mapping/mapped.tmp) | 
			# >(samtools view -@ $((${threads}/2)) -c > ${output_dir}/mapping/all.tmp) | # No need to do this now that read counting is a separate step
			samtools view -@ ${THREADS} -O bam | samtools sort -@ ${THREADS} | \
			samtools depth -aa - > ${output_dir}/coverage/${raw_read_name_base}_to_${bin_name_base}.tsv

		# Extract mapping stats
		mapped_reads=$(cat ${output_dir}/mapping/mapped.tmp)
		(>&2 echo "${mapped_reads} mapped reads")

		## Calculate additional stats - doesn't work in Bash alone because of rounding small numbers
		#percent_recruited_reads=$((${mapped_reads}/${total_reads}*100))
		#average_coverage=$((${mapped_reads}/${genome_length_nt}))
		#average_coverage_per_million_reads=$((${average_coverage}/${total_reads}*1000000))

		# Add to TSV file
		printf "${bin_name_base}\t${raw_read_name_base}\t${mapped_reads}\t${genome_length_nt}\n" >> ${bin_mapping_summary_filename}

		# Clean up
		rm ${output_dir}/mapping/mapped.tmp

		# Increase the iterator for tracking the number of processed genome-metagenome pairs
		iteration_tmp=$((${iteration}+1))
		iteration=${iteration_tmp}

	done
done

# Cleanup
rm -r ${output_dir}/mapping

(>&2 echo "[ $(date -u) ]: ${0##*/}: Finished.")

