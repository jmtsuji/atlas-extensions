refined_bin_dir=
raw_read_dir=
output_dir=


bin_paths=($(find -L ${refined_bin_dir} -iname "*.fa")) # Careful!

raw_read_paths=($(find -L ${raw_read_dir} -iname "*R1.fastq.gz")) # Grab the R1 only

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
outfile=${output_dir}/${bin_name_base}_to_${raw_read_name_base}.sam
logfile=${output_dir}/${bin_name_base}_to_${raw_read_name_base}_contig_coverage_stats.log

# Hard-code memory for now
MEMORY=50

echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: read mapping using ${#sample_filepaths_individual[@]} identified raw read files."
bbwrap.sh nodisk=t ref=${contigs} in1=${R1},${se} in2=${R2},null perfectmode=t trimreaddescriptions=t outm=${outfile} threads=${threads} pairlen=1000 pairedonly=t mdtag=t xstag=fs nmtag=t sam=1.3 local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${MEMORY}G 2> ${logfile}



done
done
