!#/usr/bin/env bash
set -euo pipefail

# Running anvi'o from atlas-coassembly output
# Jackson M. Tsuji, Neufeld lab (Apr. 10, 2018)_
# All based off the Anvi'o metagenomics tutorial at http://merenlab.org/2016/06/22/anvio-tutorial-v2/ (etc.)
# Recommended to start from within a conda env
# Still rough code!!

# User-defined variables
# TODO - should be defined as input to script
atlas_dir="/Hippodrome/jmtsuji/180123_ELA111314_atlas_r1.0.22_plus_full" # with "coassembly" dir inside
coassembly_sample_ID="CA-L227-2014"
output_dir="${atlas_dir}/post_analysis/04_anvio"
threads=12

# Hard-coded variables
#cogs_data_dir="/Hippodrome/anvio/databases/COG"
#cogs_data_dir="/Winnebago/jmtsuji/miniconda2/envs/anvio4/lib/python3.6/site-packages/anvio/data/misc/COG"

# Make needed output directories
mkdir -p ${output_dir}/01a_import_prokka ${output_dir}/01b_import_atlas_table ${output_dir}/02_multi_mapping

## Test if COGs database exists
# TODO - finish
#if [ ! -d ${cogs_data_dir}/[database??] ]; then
#	echo "ERROR: could not find COGs database at '${cogs_data_dir}/[database??]'. Exiting..."
#	exit 1
#fi

#### 1a. Import prokka info
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Exporting prokka functional gene annotations"
cd ${output_dir}/01a_import_prokka

if [ ! -f gff_parser.py ]; then
	"[$(date '+%y%m%d %H:%M:%S %Z')]: Installing gff_parser.py because it is not already present in the run folder"
	wget https://raw.githubusercontent.com/karkman/gff_parser/master/gff_parser.py -O gff_parser.py
	# TODO - make silent
fi

# TODO - check if already installed
pip install gffutils

python gff_parser.py ${atlas_dir}/coassembly/${coassembly_sample_ID}/annotation/prokka/${coassembly_sample_ID}.gff \
				--gene-calls ${coassembly_sample_ID}_gene_calls.txt \
				--annotation ${coassembly_sample_ID}_gene_annot.txt 2>&1 | tee gff_parser.log


#### 1b. Import ATLAS table info
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Exporting information from the ATLAS annotations table"
cd ${output_dir}/01b_import_atlas_table

if [ ! -f parse_atlas_table_for_anvio.R ]; then
	# TODO - do this with proper version control
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Installing parse_atlas_table_for_anvio.R because it is not already present in the run folder"
	git clone https://github.com/jmtsuji/atlas-extensions.git
	mv atlas-extensions/parse_atlas_table_for_anvio.R .
	rm -rf atlas-extensions
	chmod 755 parse_atlas_table_for_anvio.R
fi

./parse_atlas_table_for_anvio.R -a ${atlas_dir}/coassembly/${coassembly_sample_ID}/${coassembly_sample_ID}_annotations_multi_mapped.txt \
				-t ${coassembly_sample_ID}_gene_taxonomy.tsv -c ${coassembly_sample_ID}_binning_results.tsv \
				-b ${coassembly_sample_ID}_bins_info.tsv 2>&1 | tee parse_atlas_table_for_anvio.log


#### 2. Generate contigs database
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Generating the contigs database"
cd ${output_dir}
anvi-gen-contigs-database -f ${atlas_dir}/coassembly/${coassembly_sample_ID}_contigs.fasta \
				-o ${coassembly_sample_ID}_contigs.db -n ${coassembly_sample_ID}_contigs_db \
				--external-gene-calls ${output_dir}/01a_import_prokka/${coassembly_sample_ID}_gene_calls.txt

## Example
# gene_callers_id	contig	start	stop	direction	partial	source	version
# 1	contig_01	1113	1677	f	0	program	v1.0
## "The statement above means that the index of the first nucleotide in any contig should be 0. In other words, we start counting from 0, not from 1."


#### 3. Annotate with single-copy marker genes
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Annotating with single-copy marker genes"
cd ${output_dir}
anvi-run-hmms -c ${coassembly_sample_ID}_contigs.db --num-threads ${threads}
# anvi-run-hmms -c ${coassembly_sample_ID}_contigs.db --num-threads ${threads} --hmm-profile-dir [your_dir_for_custom_hmms]
# anvi-display-contigs-stats ${coassembly_sample_ID}_contigs.db # Will this work?


#### 4. Add functional info
# TODO - find a way to test whether or not setup is needed. Assumes already set up for now.
# anvi-setup-ncbi-cogs --num-threads ${threads} --cog-data-dir ${cogs_data_dir} --cog-data-source ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Annotating with COGs"
cd ${output_dir}
anvi-run-ncbi-cogs --cog-data-dir ${cogs_data_dir} --num-threads ${threads}


#### 5. Import functional and taxonomic info
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing functional annotations from prokka"
cd ${output_dir}
anvi-import-functions -c ${coassembly_sample_ID}_contigs.db -i ${output_dir}/01a_import_prokka/${coassembly_sample_ID}_gene_annot.txt

## Example table
# gene_callers_id	source	accession	function	e_value
# 1	Pfam	PF01132	Elongation factor P (EF-P) OB domain	4e-23

# TODO - format ATLAS annotations table into the input matrix. Need to parse Greengenes.
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing taxonomic gene classifications from ATLAS"
anvi-import-taxonomy -c ${coassembly_sample_ID}_contigs.db -i 01b_import_atlas_table/${coassembly_sample_ID}_gene_taxonomy.tsv -p default_matrix

## Example table
# gene_callers_id	t_phylum	t_class	t_order	t_family	t_genus	t_species
# 1	 	 	 	 	 	Bacteroides fragilis
# 2	 	 	 	 	 	Bacteroides fragilis


#### 6. Make relative abundance profiles
# Get samples
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Mapping relative abundance profiles"
cd ${output_dir}/02_mapping
sample_names=($(find ${atlas_dir}/coassembly/${coassembly_sample_ID}/multi_mapping -name "*.bam" -type f))

for sample in ${sample_names[@]}; do
	sample_name=${sample##*/}
	sample_name=${sample_name%.*}
	
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Mapping ${sample_name}"
	anvi-profile -i ${sample} -c ${coassembly_sample_ID}_contigs.db --output-dir ${sample_name} --sample-name ${sample_name}
done

echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Merging information from mapped samples"
cd ${output_dir}
anvi-merge ${output_dir}/02_mapping/*/PROFILE.db -o ${coassembly_sample_ID}_samples_merged -c ${coassembly_sample_ID}_contigs.db --skip-concoct-binning

#### 7. Import my bins
# TODO - get my bin info into the proper format - 2 tables needed.

echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing bin info from ATLAS"
cd ${output_dir}
anvi-import-collection 01b_import_atlas_table/${coassembly_sample_ID}_binning_results.tsv -p ${coassembly_sample_ID}_samples_merged/PROFILE.db \
				-c ${coassembly_sample_ID}_contigs.db --source "metabat2" --contigs-mode \
				--bins-info 01b_import_atlas_table/${coassembly_sample_ID}_bins_info.tsv

## Example table: external_binning_of_contigs.txt
# 204_10M_MERGED.PERFECT.gz.keep_contig_878	Bin_2
# 204_10M_MERGED.PERFECT.gz.keep_contig_6515	Bin_3
# 204_10M_MERGED.PERFECT.gz.keep_contig_1720	Bin_4
# contig_1	Bin_4

## Example info table: example_bins_info_file.txt
# Bin_0	UNKNOWN_SOURCE	#ABCDEF
# Bin_1	UNKNOWN_SOURCE	#FF2244
# Bin_2	UNKNOWN_SOURCE	#2244FF
# Bin_3	UNKNOWN_SOURCE	#22FF44
# Bin_4	UNKNOWN_SOURCE	#44FFFF

#### 8. Visualize - run OUTSIDE shell script.
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Pipeline finished."
printf "\nTo visualize, please run:\n"
printf "cd ${output_dir}\n"
printf "anvi-interactive -p ${coassembly_sample_ID}_samples_merged/PROFILE.db -c ${coassembly_sample_ID}_contigs.db -C metabat2 --server-only -P 8080\n\n"

