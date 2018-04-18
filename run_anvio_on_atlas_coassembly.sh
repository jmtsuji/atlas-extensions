#!/usr/bin/env bash
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
threads=6

# Hard-coded variables
#cogs_data_dir="/Hippodrome/anvio/databases/COG"
#cogs_data_dir="/Winnebago/jmtsuji/miniconda2/envs/anvio4/lib/python3.6/site-packages/anvio/data/misc/COG"

function export_prokka_info {
	cd ${output_dir}/01a_import_prokka

	if [ ! -f gff_parser.py ]; then
		
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Installing gff_parser.py because it is not already present in the run folder."
		wget --quiet -O gff_parser.py https://raw.githubusercontent.com/karkman/gff_parser/master/gff_parser.py
		# Consider removing '--quiet' flag from wget in case helpful for logging
		
	fi
	
	# Install gffutils if not installed
	if pip list 2>/dev/null | grep -q gffutils; then

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Python package 'gffutils' is already installed; will not re-install."	
		
	else
		
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Installing python package 'gffutils' via pip because it is not already installed."
		pip install --quiet gffutils

	fi
		
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Exporting prokka functional gene annotations"
	python gff_parser.py ${atlas_dir}/coassembly/${coassembly_sample_ID}/annotation/prokka/${coassembly_sample_ID}.gff \
					--gene-calls ${coassembly_sample_ID}_gene_calls.txt --annotation ${coassembly_sample_ID}_gene_annot.txt
	# Script gives no log info.

}

function export_atlas_info {

	cd ${output_dir}/01b_import_atlas_table

	if [ ! -f parse_atlas_table_for_anvio.R ]; then

		# TODO - do this with proper version control
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Installing parse_atlas_table_for_anvio.R because it is not already present in the run folder"
		git clone --quiet https://github.com/jmtsuji/atlas-extensions.git
		mv atlas-extensions/parse_atlas_table_for_anvio.R .
		rm -rf atlas-extensions
		chmod 755 parse_atlas_table_for_anvio.R
		
	else

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: 'parse_atlas_table_for_anvio.R' is already present in the run folder; will not re-install."
		
	fi

	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Exporting information from the ATLAS annotations table"
	./parse_atlas_table_for_anvio.R -a ${atlas_dir}/coassembly/${coassembly_sample_ID}/${coassembly_sample_ID}_annotations_multi_mapped.txt \
					-t ${coassembly_sample_ID}_gene_taxonomy.tsv -c ${coassembly_sample_ID}_binning_results.tsv \
					-b ${coassembly_sample_ID}_bins_info.tsv 2>&1 | tee parse_atlas_table_for_anvio.log

	# TODO - deal with the differing length of the gene calls from the gff file versus the gene taxonomy

}

function generate_contig_database {

	#### 2. Generate contigs database
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Generating the contigs database"
	cd ${output_dir}
	anvi-gen-contigs-database -f ${atlas_dir}/coassembly/${coassembly_sample_ID}/${coassembly_sample_ID}_contigs.fasta \
					-o ${coassembly_sample_ID}_contigs.db -n ${coassembly_sample_ID}_contigs_db \
					--external-gene-calls ${output_dir}/01a_import_prokka/${coassembly_sample_ID}_gene_calls.txt \
					--ignore-internal-stop-codons --split-length -1 2>&1 | tee misc_logs/anvi-gen-contigs-database.log

	## Example table
	# gene_callers_id	contig	start	stop	direction	partial	source	version
	# 1	contig_01	1113	1677	f	0	program	v1.0
	## "The statement above means that the index of the first nucleotide in any contig should be 0. In other words, we start counting from 0, not from 1."

}

function add_hmm_annotations {

	#### Annotate with single-copy marker genes
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Annotating with single-copy marker genes"
	cd ${output_dir}
	anvi-run-hmms -c ${coassembly_sample_ID}_contigs.db --num-threads ${threads} \
					2>&1 | tee misc_logs/anvi-run-hmms.log
	# anvi-run-hmms -c ${coassembly_sample_ID}_contigs.db --num-threads ${threads} --hmm-profile-dir [your_dir_for_custom_hmms]
	# anvi-display-contigs-stats ${coassembly_sample_ID}_contigs.db # Will this work?

}

function add_cog_annotations {

	# TODO - find a way to test whether or not setup is needed. Assumes already set up for now.
	
	# To set up:
	# anvi-setup-ncbi-cogs --num-threads ${threads} --just-do-it 2>&1 | tee misc_logs/anvi-setup-ncbi-cogs.log # --cog-data-dir ${cogs_data_dir}
	
	## To test if custom COGs database exists:
	# TODO - finish or delete
	#if [ ! -d ${cogs_data_dir}/[database??] ]; then
	#	echo "ERROR: could not find COGs database at '${cogs_data_dir}/[database??]'. Exiting..."
	#	exit 1
	#fi
	
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Annotating with COGs"
	cd ${output_dir}
	mkdir tmp
	anvi-run-ncbi-cogs --num-threads ${threads} -c ${coassembly_sample_ID}_contigs.db \
					--temporary-dir-path tmp 2>&1 | tee misc_logs/anvi-run-ncbi-cogs.log # --cog-data-dir ${cogs_data_dir}
	rm -rf tmp

}

function import_atlas_annotations {

	#### Import functional and taxonomic info
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing functional annotations from prokka"
	cd ${output_dir}
	anvi-import-functions -c ${coassembly_sample_ID}_contigs.db \
					-i ${output_dir}/01a_import_prokka/${coassembly_sample_ID}_gene_annot.txt \
					2>&1 | tee misc_logs/anvi-import-functions.log

	## Example table
	# gene_callers_id	source	accession	function	e_value
	# 1	Pfam	PF01132	Elongation factor P (EF-P) OB domain	4e-23

	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing taxonomic gene classifications from ATLAS"
	anvi-import-taxonomy -c ${coassembly_sample_ID}_contigs.db \
					-i 01b_import_atlas_table/${coassembly_sample_ID}_gene_taxonomy.tsv \
					-p default_matrix 2>&1 | tee misc_logs/anvi-import-taxonomy.log

	## Example table
	# gene_callers_id	t_phylum	t_class	t_order	t_family	t_genus	t_species
	# 1	 	 	 	 	 	Bacteroides fragilis
	# 2	 	 	 	 	 	Bacteroides fragilis

}

function make_read_mapping_profiles {

	# Get samples
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Mapping relative abundance profiles"
	cd ${output_dir}/02_multi_mapping
	sample_names=($(find ${atlas_dir}/coassembly/${coassembly_sample_ID}/multi_mapping -name "*.bam" -type f))

	for sample in ${sample_names[@]}; do
		sample_name=${sample##*/}
		sample_name=${sample_name%.*}
		sample_name_simple=$(sed s/'-'/'_'/g <<<${sample_name}) # Get rid of dashes
		
		# Index BAM
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: sorting and indexing BAM"
		#anvi-init-bam -o ${sample_name_simple}.bam ${sample}
		samtools sort -o ${sample_name_simple}.bam -@ ${threads} -m 4G ${sample}
		samtools index -b -@ ${threads} ${sample_name_simple}.bam
		
		# Generate profile
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: creating mapping profile"
		anvi-profile -i ${sample_name_simple}.bam -c ${output_dir}/${coassembly_sample_ID}_contigs.db \
						--output-dir ${sample_name_simple} --sample-name ${sample_name_simple} -T ${threads} \
						--min-contig-length 1000 2>&1 | tee logs/anvi-profile_${sample_name_simple}.log
		
		# TODO - uncomment this once the code seems to be working reliably.
		## Remove indexed bam
		# rm ${sample_name_simple}.bam ${sample_name_simple}.bam.bai
		
	done

}

function merge_read_mapping_profiles {

	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Merging information from mapped samples into contig database"
	cd ${output_dir}
	anvi-merge ${output_dir}/02_multi_mapping/*/PROFILE.db -o ${coassembly_sample_ID}_samples_merged \
					-c ${coassembly_sample_ID}_contigs.db --skip-concoct-binning --S metabat2 \
					2>&1 | tee misc_logs/anvi-merge.log

}

function import_custom_bins {

	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing bin info from ATLAS"
	cd ${output_dir}
	anvi-import-collection 01b_import_atlas_table/${coassembly_sample_ID}_binning_results.tsv \
					-p ${coassembly_sample_ID}_samples_merged/PROFILE.db \
					-c ${coassembly_sample_ID}_contigs.db -C "metabat2" --contigs-mode \
					--bins-info 01b_import_atlas_table/${coassembly_sample_ID}_bins_info.tsv \
					2>&1 | tee misc_logs/anvi-import-collection.log

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

}

function main {

# Make needed output directories
mkdir -p ${output_dir}/01a_import_prokka ${output_dir}/01b_import_atlas_table ${output_dir}/02_multi_mapping/logs ${output_dir}/misc_logs

export_prokka_info
export_atlas_info

generate_contig_database

add_hmm_annotations
add_cog_annotations
import_atlas_annotations

make_read_mapping_profiles
merge_read_mapping_profiles

import_custom_bins

# Exit and report anvi-interative instructions to user
echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Pipeline finished."
printf "\nTo visualize, please run:\n"
printf "cd ${output_dir}\n"
printf "anvi-interactive -p ${coassembly_sample_ID}_samples_merged/PROFILE.db -c ${coassembly_sample_ID}_contigs.db -C metabat2 --server-only -P 8080\n\n"

}

main

