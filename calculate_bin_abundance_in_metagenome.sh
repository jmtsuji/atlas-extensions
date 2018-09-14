#!/usr/bin/env bash
set -euo pipefail

# If no input is provided, provide help and exit
if [ $# == 0 ]
    then
    printf "$(basename $0): Wrapper to calculate the read recruitment of genome bins to metagenomes.\n"
    printf "Contact Jackson M. Tsuji (jackson.tsuji@uwaterloo.ca) for bug reports or feature requests.\n\n"
    printf "Installation: you can install all dependencies via conda and then add the scripts from the Github repo with:\n"
    printf "              conda create -y -n calc_bin_abundance -c conda-forge -c bioconda samtools bbmap\n"
    printf "              source activate calc_bin_abundance\n\n"
    printf "Usage: $(basename $0) refined_bin_dir raw_read_dir output_dir threads memory 2>&1 | tee $(basename $0 .sh).log\n\n"
    printf "Usage details:\n"
    printf "   1. refined_bin_dir: Directory containing genomes/bins. MUST have extension *.fa!! Symlinks are okay.\n"
    printf "   2. raw_read_dir: Directory containing QC'ed, unassembled metagenomic reads with naming structure matching that of the ATLAS pipeline.\n"
    printf "                    *QC_R1.fastq.gz, *QC_R2.fastq.gz, and *QC_se.fastq.gz all needed. Symlinks are okay.\n"
    printf "   3. output_dir: Directory where you want the results to be output. Anything already there might be overwritten.\n"
    printf "   4. threads: number of threads to run.\n"
    printf "   5. memory: in gigabytes.\n\n"
    exit 1
fi

# Set user variables
refined_bin_dir=$1 # bins MUST have extension *.fa
raw_read_dir=$2 # reads MUST be output of ATLAS QC and match naming structure. R1, R2, se reads all need to be there, for all samples.
output_dir=$3
THREADS=$4
MEMORY=$5 # in Gigabytes

(>&2 echo "[ $(date -u) ]: Running ${0##*/}")
(>&2 echo "[ $(date -u) ]: refined_bin_dir: ${refined_bin_dir}")
(>&2 echo "[ $(date -u) ]: raw_read_dir: ${raw_read_dir}")
(>&2 echo "[ $(date -u) ]: output_dir: ${output_dir}")
(>&2 echo "[ $(date -u) ]: bin_mapping_summary_filename: bin_mapping_stats.tsv (in output_dir)")
(>&2 echo "[ $(date -u) ]: threads: ${THREADS}")
(>&2 echo "[ $(date -u) ]: memory: ${MEMORY} GB")

# Create folder structure in output dir
mkdir -p ${output_dir}/mapping ${output_dir}/logs
bin_mapping_summary_filename="${output_dir}/bin_mapping_stats.tsv"

# Find bin and raw read files
bin_paths=($(find -L ${refined_bin_dir} -iname "*.fa")) # Careful - make sure the files match *.fa!
raw_read_paths=($(find -L ${raw_read_dir} -iname "*R1.fastq.gz")) # Grab the R1 only
iterations=$((${#bin_paths[@]}*${#raw_read_paths[@]}))

(>&2 echo "[ $(date -u) ]: found ${#bin_paths[@]} genome bins with extension *.fa in the refined_bin_dir")
(>&2 echo "[ $(date -u) ]: found ${#raw_read_paths[@]} set of unassembled metagenomic reads in the raw_read_dir")
(>&2 echo "[ $(date -u) ]: this means found ${#bin_paths[@]} * ${#raw_read_paths[@]} = ${iterations} combinations of mapping will be run")

# Initialize output read stats file
(>&2 echo "[ $(date -u) ]: Initiatilizing '${bin_mapping_summary_filename##*/}'")
printf "mapping_filename\tgenome\tmetagenome\tmapped_reads\ttotal_reads\tgenome_length_nt\n" > ${bin_mapping_summary_filename}

# Run all against all in a nested 'for' loop
for bin_path in ${bin_paths[@]}; do

	# Get total nucleotide length in bin for later
	genome_char_count=$(grep -v "^>" ${bin_path} | wc)
	genome_new_lines=$(echo ${genome_char_count} | cut -d ' ' -f 1)
	genome_total_chars=$(echo ${genome_char_count} | cut -d ' ' -f 3)
	genome_length_nt=$((${genome_total_chars}-${genome_new_lines}))

	(>&2 echo "[ $(date -u) ]: Starting on bin '${bin_path##*/}' with length ${genome_length_nt} nt")

	for raw_read_path in ${raw_read_paths[@]}; do

		# Get base names of bins and reads
		bin_name_base=${bin_path%.*}
		bin_name_base=${bin_name_base##*/}
		raw_read_name_base=${raw_read_path%*_QC_R1.fastq.gz}
		raw_read_name_base=${raw_read_name_base##*/}

		# Assign generic variables for read mapping
		R1=${raw_read_dir}/${raw_read_name_base}_QC_R1.fastq.gz
		R2=${raw_read_dir}/${raw_read_name_base}_QC_R2.fastq.gz
		se=${raw_read_dir}/${raw_read_name_base}_QC_se.fastq.gz
		contigs=${bin_path}
		read_mapping_file=${output_dir}/mapping/${bin_name_base}_to_${raw_read_name_base}.sam
		logfile=${output_dir}/logs/${bin_name_base}_to_${raw_read_name_base}_contig_coverage_stats.log

		# Read map AND pipe directly to stats (to avoid excessive input/output, which is rough on hard drives)
		(>&2 echo "[ $(date -u) ]: mapping '${raw_read_name_base}*fastq.gz' to '${bin_name_base}' (log: logs/${bin_name_base}_to_${raw_read_name_base}_contig_coverage_stats.log)")
		bbwrap.sh nodisk=t ref=${contigs} in1=${R1},${se} in2=${R2},null perfectmode=t trimreaddescriptions=t \
			out=stdout threads=${THREADS} pairlen=1000 pairedonly=t mdtag=t xstag=fs nmtag=t sam=1.3 \
			local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${MEMORY}G 2> ${logfile} | \
			tee >(samtools view -@ $((${THREADS}/2)) -c -F 4 > ${output_dir}/mapping/mapped.tmp) | 
			samtools view -@ $((${THREADS}/2)) -c > ${output_dir}/mapping/all.tmp

		# (>&2 echo "[ $(date -u) ]: ${raw_read_name_base} to ${bin_name_base}: Calculating stats")

		# Extract mapping stats
		mapped_reads=$(cat ${output_dir}/mapping/mapped.tmp)
		total_reads=$(cat ${output_dir}/mapping/all.tmp)
		(>&2 echo "[ $(date -u) ]: ${raw_read_name_base} to ${bin_name_base}: ${mapped_reads} mapped; ${total_reads} total. Adding to TSV.")

		## Calculate additional stats
		#percent_recruited_reads=$((${mapped_reads}/${total_reads}*100))
		#average_coverage=$((${mapped_reads}/${genome_length_nt}))
		#average_coverage_per_million_reads=$((${average_coverage}/${total_reads}*1000000))
		
		# Add to TSV file
		printf "${read_mapping_file##*/}\t${bin_name_base}\t${raw_read_name_base}\t${mapped_reads}\t${total_reads}\t${genome_length_nt}\n" >> ${bin_mapping_summary_filename}

		# Clean up
		rm ${output_dir}/mapping/mapped.tmp ${output_dir}/mapping/all.tmp

	done
done

(>&2 echo "[ $(date -u) ]: ${0##*/}: Finished.")

