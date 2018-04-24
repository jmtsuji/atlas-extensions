#!/usr/bin/env bash
set -euo pipefail
script_version="1.0.22-coassembly-r4-dev" # to match ATLAS version it is designed to work with

# atlas-to-anvi
# Copyright Jackson M. Tsuji, 2018
# Neufeld lab, University of Waterloo, Canada
# Created Apr. 10, 2018
# Description: Imports ATLAS-coassembly output into anvio for manual bin refinement, for a specific coassembly sample.
# REQUIRES that ATLAS coassembly (via atlas-extensions version 1.0.22-coassembly-r3) has been run and that the desired sample is in the 'coassembly' directory within the ATLAS run folder.
# All based off the Anvi'o metagenomics tutorial at http://merenlab.org/2016/06/22/anvio-tutorial-v2/ (etc.)

# # Recommended to start from within a conda env, e.g.,
# conda create -y -n anvio4 python=3
# conda install -y -n anvio4 -c bioconda -c conda-forge anvio=4
# source activate anvio4

# If no input is provided, provide help and exit
if [ $# == 0 ]
    then
    printf "$(basename $0): Imports ATLAS-coassembly output into anvio for manual bin refinement. Early development version.\n"
    printf "Version: ${script_version}\n"
    printf "Contact Jackson M. Tsuji (jackson.tsuji@uwaterloo.ca) for bug reports or feature requests.\n\n"
    printf "Usage: $(basename $0) atlas_dir coassembly_sample_ID output_dir threads 2>&1 | tee $(basename $0 .sh).log\n\n"
    printf "Usage details:\n"
    printf "1. atlas_dir: Path to the base directory where ATLAS files were output.\n"
    printf "2. coassembly_sample_ID: Exact name of the coassembly that you desire to visualize bins from. MUST be a sample within the 'coassembly' directory in the 'atlas_dir' folder.\n"
    printf "3. output_dir: Path to save the anvio database and associated files to.\n"
    printf "4. threads: number of threads to run.\n\n"
    printf "Additional usage notes:\n"
    printf "* REQUIREMENTS: supports coassembly output from atlas-extensions version '1.0.22-coassembly-r3', although you do not need to have the requirements for the ATLAS coassembly extension to run this script. Must have anvio4 installed (ideally via conda). Must have the Anvi'o COG database set up via 'anvi-setup-ncbi-cogs'. Also, need R and an internet connection.\n"
    exit 1
fi

# Example
# atlas-to-anvi.sh "/Hippodrome/jmtsuji/180123_ELA111314_atlas_r1.0.22_plus_full" "CA-L227-2013" "/Hippodrome/jmtsuji/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/04_anvio/CA-L227-2013" 6 2>&1 | tee atlas-to-anvi.log

# Set variables from user input:
atlas_dir=$1
coassembly_sample_ID=$2
output_dir=$3
threads=$4

function test_inputs {

	if [ ! -d ${atlas_dir} ]; then
	
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: atlas_dir directory '${atlas_dir}' does not exist. Is your atlas_dir correct? Exiting..."
		exit 1
	
	fi
	
	
	if [ ! -d ${atlas_dir}/coassembly ]; then
	
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: directory '${atlas_dir}/coassembly' does not exist. Have you run atlas-coassembly? Exiting..."
		exit 1
	
	fi


	if [ ! -d ${atlas_dir}/coassembly/${coassembly_sample_ID} ]; then
	
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: directory '${atlas_dir}/coassembly/${coassembly_sample_ID}' does not exist. Is your coassembly_sample_ID '${coassembly_sample_ID}' correct? Exiting..."
		exit 1
	
	fi
	
	
	if [ -d ${output_dir} ]; then
	
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: directory '${output_dir}' already exists. Please delete before starting this script (if you don't care about the data in that folder), or choose a different output directory. Exiting..."
		exit 1
	
	fi

}

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
		
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Exporting prokka gene annotations"
	python gff_parser.py ${atlas_dir}/coassembly/${coassembly_sample_ID}/annotation/prokka/${coassembly_sample_ID}.gff \
					--gene-calls ${coassembly_sample_ID}_gene_calls.txt --annotation ${coassembly_sample_ID}_gene_annot.txt
	# Script gives no log info.

}

## TODO - add a different sanity check. Function in its current form does not work if genes without taxonomic assignment are at the end of the file...
function match_atlas_table_to_prokka_info {
	# TODO - move variables to local (within function) instead of global
	
	cd ${output_dir}/01b_import_atlas_table
	
	# Get the final gene entry number in the prokka file
	last_prokka_gene_ID=$(tail -n 1 ${output_dir}/01a_import_prokka/${coassembly_sample_ID}_gene_calls.txt | cut -d $'\t' -f 1)
	
	# Get the final gene entry number of the altas taxonomy file
	last_atlas_gene_ID=$(tail -n 1 ${output_dir}/01b_import_atlas_table/${coassembly_sample_ID}_gene_taxonomy.tsv | cut -d $'\t' -f 1)

	if [ ${last_prokka_gene_ID} -eq ${last_atlas_gene_ID} ]; then
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Prokka gene calls and ATLAS gene taxonomy annotations match in total numbers - good."
		
	elif [ ${last_prokka_gene_ID} -lt ${last_atlas_gene_ID} ]; then
		
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: WARNING: prokka gene calls and ATLAS gene taxonomy annotations do NOT match in total numbers. Final prokka gene call ID is '${last_prokka_gene_ID}', but final ATLAS taxonomy gene annotation ID is '${last_atlas_gene_ID}'. This has been observed before possibly as an error in prokka while generating the GFF file (final entries are cut out?)."

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: To correct error, deleting gene IDs in ATLAS gene taxonomy file from after ${last_prokka_gene_ID} until ${last_atlas_gene_ID}. Keeping old (complete) gene taxonomy file as '${coassembly_sample_ID}_gene_taxonomy_ORIGINAL.tsv' but saving the new (truncated) gene taxonomy file over the original at '${coassembly_sample_ID}_gene_taxonomy.tsv'."
		
		# Find the line number in the altas gene taxonomy file where the 'last_prokka_gene_ID' is
		atlas_tax_file="${output_dir}/01b_import_atlas_table/${coassembly_sample_ID}_gene_taxonomy.tsv"
		
		# First, do a sanity check that there is only one match to the prokka gene caller ID
		num_hits=$(grep -c ^${last_prokka_gene_ID} ${atlas_tax_file})
		if [ ${num_hits} -eq 1 ]; then
			# Then get the line number if all looks okay.
			matching_atlas_tax_line=$(grep -n ^${last_prokka_gene_ID} ${atlas_tax_file} | cut -d ":" -f 1)
		else
			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: Found duplicated entries in the atlas taxonomy table '${atlas_tax_file}'. ${num_hits} rows begin with '${last_prokka_gene_ID}' (should only be one row). Exiting..."
			exit 1
		fi
		
		# Keep backup of original taxonomy file but make a new truncated one for the remainder of the pipeline.
		mv ${coassembly_sample_ID}_gene_taxonomy.tsv ${coassembly_sample_ID}_gene_taxonomy_ORIGINAL.tsv
		head -n ${matching_atlas_tax_line} ${coassembly_sample_ID}_gene_taxonomy_ORIGINAL.tsv > ${coassembly_sample_ID}_gene_taxonomy.tsv
		
	elif [ ${last_prokka_gene_ID} -gt ${last_atlas_gene_ID} ]; then
	
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: prokka gene calls and ATLAS gene taxonomy annotations do NOT match in total numbers. Final prokka gene call ID is '${last_prokka_gene_ID}', but final ATLAS taxonomy gene annotation ID is '${last_atlas_gene_ID}'. The case where the ATLAS annotations are truncated compared to the prokka annotations has not been seen before. Something could be wrong with your ALTAS inputs, or this could be an unknown bug. Exiting..."
		exit 1
	
	fi
	
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
	./parse_atlas_table_for_anvio.R -a ${atlas_dir}/coassembly/${coassembly_sample_ID}/${coassembly_sample_ID}_annotations_multi_mapping.txt \
					-t ${coassembly_sample_ID}_gene_taxonomy.tsv -c ${coassembly_sample_ID}_binning_results.tsv \
					-b ${coassembly_sample_ID}_bins_info.tsv 2>&1 -@ ${threads} | tee parse_atlas_table_for_anvio.log

	# TODO - add a different sanity check. Function in its current form does not work if genes without taxonomic assignment are at the end of the file...
	## Deal with the differing length of the gene calls from the gff file versus the gene taxonomy
	# match_atlas_table_to_prokka_info

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
	cd ${output_dir}

	# TODO - get past the error within anvi'o where this seems incompatible with COG annotations. Skipping for now.
	#### Import functional info
	#echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing functional annotations from prokka"
	#anvi-import-functions -c ${coassembly_sample_ID}_contigs.db \
	#				-i ${output_dir}/01a_import_prokka/${coassembly_sample_ID}_gene_annot.txt \
	#				2>&1 | tee misc_logs/anvi-import-functions.log

	## Example table
	# gene_callers_id	source	accession	function	e_value
	# 1	Pfam	PF01132	Elongation factor P (EF-P) OB domain	4e-23

	
	#### Import taxonomy info
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
		#anvi-init-bam -o ${sample_name_simple}.bam ${sample} # TODO - delete this step. Can use samtools directly instead for better efficiency.
		#samtools sort -o ${sample_name_simple}.bam -@ ${threads} -m 4G ${sample} # TODO - sort already done in ATLAS. Delete this step.
		
		cp ${sample} ${sample_name_simple}.bam
		samtools index -b -@ ${threads} ${sample_name_simple}.bam
		# TODO - 'cp' step above: consider a more efficient way to do this long-term. I would just run the index step on the original ${sample}, but it might be in a write-protected directory. However, it seems inefficient to copy such large files (this is also hard on the disk).
		
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
					-c ${coassembly_sample_ID}_contigs.db --skip-concoct-binning -S metabat2 \
					2>&1 | tee misc_logs/anvi-merge.log
	# Consider '--enforce-hierarchical-clustering' to cluster even with > 25,000 contigs. But could take a long time...

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

function summarize {

	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Summarizing anvio binning statistics"
	cd ${output_dir}
	anvi-summarize -p ${coassembly_sample_ID}_samples_merged/PROFILE.db \
					-c ${coassembly_sample_ID}_contigs.db -C metabat2 \
					-o ${coassembly_sample_ID}_summary --taxonomic-level t_genus \
					--init-gene-coverages 2>&1 | tee misc_logs/anvi-summarize.log
	# Note: '--init-gene-coverages' takes a long time to run.
	# Note: consider '--quick-summary' for faster run with with more minimal output.

}

function main {

	echo "Running $(basename $0), version ${script_version}, on $(date)."
	echo ""

	# Get date and time of start
	start_time=$(date)
	
	# Report settings
	echo "Run settings:"
	echo "atlas_dir: ${atlas_dir}"
	echo "coassembly_sample_ID: ${coassembly_sample_ID}"
	echo "output_dir: ${output_dir}"
	echo "threads: ${threads}"
	echo ""

	test_inputs
	
	# Make needed output directories
	mkdir -p ${output_dir}/01a_import_prokka ${output_dir}/01b_import_atlas_table ${output_dir}/02_multi_mapping/logs ${output_dir}/misc_logs

	# Get relevant annotations from the ATLAS run
	export_prokka_info
	export_atlas_info

	# Create the database
	generate_contig_database

	# Annotate the database
	add_hmm_annotations
	add_cog_annotations
	import_atlas_annotations

	# Add read mapping information to the database
	make_read_mapping_profiles
	merge_read_mapping_profiles

	# Import bins
	import_custom_bins
	
	# Generate summary
	summarize
	
	end_time=$(date)

	# Exit and report anvi-interative instructions to user
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: $(basename $0): complete. Started at ${start_time} and finished at ${end_time}."
	printf "\nTo visualize, please run:\n"
	printf "cd ${output_dir}\n"
	printf "anvi-interactive -p ${coassembly_sample_ID}_samples_merged/PROFILE.db -c ${coassembly_sample_ID}_contigs.db -C metabat2 --server-only -P 8080\n"
	printf "bin_name=[name_of_interest]\n"
	printf "anvi-refine -p ${coassembly_sample_ID}_samples_merged/PROFILE.db -c ${coassembly_sample_ID}_contigs.db -C metabat2 -b \${bin_name} --taxonomic-level t_genus --title \${bin_name} -P 8080 --server-only\n\n"

}

main

