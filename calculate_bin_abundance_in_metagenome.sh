#!/usr/bin/env bash
set -euo pipefail

## To install needed dependendies:
# conda create -y -n calc_bin_abundance -c conda-forge -c bioconda samtools bbmap
# source calc_bin_abundance

# Set user variables (hard-coded for now)
# TODO - allow user to set themselves
refined_bin_dir= # bins MUST have extension *.fa
raw_read_dir= # reads MUST be output of ATLAS QC and match naming structure. R1, R2, se reads all need to be there, for all samples.
output_dir=
THREADS=40
MEMORY=50

# Create folder structure in output dir
mkdir -p ${output_dir}/mapping ${output_dir}/logs
output_filename="${output_dir}/bin_mapping_stats.tsv"

# Find bin and raw read files
bin_paths=($(find -L ${refined_bin_dir} -iname "*.fa")) # Careful - make sure the files match *.fa!
raw_read_paths=($(find -L ${raw_read_dir} -iname "*R1.fastq.gz")) # Grab the R1 only

# Initialize output read stats file
# TODO - add column for average coverage
printf 'mapping_filename\tgenome\tmetagenome\tmapped_reads\ttotal_reads\n' > ${output_filename}

# Run all against all in a nested for loop
for bin_path in ${bin_paths}; do
for raw_read_path in ${raw_read_paths}; do

# Get base names of bins and reads
bin_name_base=${bin_path%.*}
bin_name_base=${bin_name_base##*/}
raw_read_name_base=${raw_read_path%*R1.fastq.gz}
raw_read_name_base=${raw_read_name_base##*/}

# Assign generic variables for read mapping
R1=${raw_read_dir}/${raw_read_name_base}R1.fastq.gz
R2=${raw_read_dir}/${raw_read_name_base}R2.fastq.gz
se=${raw_read_dir}/${raw_read_name_base}se.fastq.gz
contigs=${bin_path}
outfile=${output_dir}/mapping/${bin_name_base}_to_${raw_read_name_base}.sam
logfile=${output_dir}/logs/${bin_name_base}_to_${raw_read_name_base}_contig_coverage_stats.log

# Read map
(>&2 echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: mapping '${raw_read_name_base}*.fastq.gz' to '${bin_name_base}'.")
bbwrap.sh nodisk=t ref=${contigs} in1=${R1},${se} in2=${R2},null perfectmode=t trimreaddescriptions=t \
	out=${outfile} threads=${THREADS} pairlen=1000 pairedonly=t mdtag=t xstag=fs nmtag=t sam=1.3 \
	local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${MEMORY}G 2> ${logfile}

(>&2 echo "[ $(date -u) ]: Calculating stats for ${mapping_filename##*/}")

# Calculate mapping stats
# TODO - add average coverage
mapped_reads=$(samtools view -@ ${THREADS} -c -F 4 ${outfile})
total_reads=$(samtools view -@ ${THREADS} -c ${outfile})
(>&2 echo "[ $(date -u) ]: ${outfile##*/}: ${mapped_reads} mapped; ${total_reads} total.")

# Add to TSV file
# TODO - add average coverage
printf "${outfile##*/}\t${bin_name_base}\t${raw_read_name_base}\t${mapped_reads}\t${total_reads}\n" >> ${output_filename}

done
done

(>&2 echo "[ $(date -u) ]: Finished.")
