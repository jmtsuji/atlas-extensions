#!/usr/bin/env bash
set -euo pipefail
script_version="1.0.22-coassembly-r4-dev2" # to match ATLAS version it is designed to work with

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
    printf "Installation: you can install all dependencies via conda and then add the scripts from the Github repo with:\n"
    printf "              conda create -y -n atlas_to_anvi -c bioconda -c conda-forge -c r anvio diamond bwa bbmap gffutils r r-plyr r-dplyr r-getopt\n"
    printf "              git clone https://github.com/jmtsuji/atlas-extensions.git\n"
    printf "              cp atlas-extensions/parse*R atlas-extensions/atlas-to-anvi.sh \${HOME}/miniconda*/envs/atlas_to_anvi/bin\n"
    printf "              rm -rf atlas-extensions\n"
    printf "              source activate atlas_to_anvi\n"
    printf "              anvi-setup-ncbi-cogs --just-do-it # run the first time you ever set up the environment"
    printf "              # Then, you're ready to go!.\n\n"
    printf "Usage: $(basename $0) run_mode atlas_dir assembly_sample_ID output_dir threads mapping_guide_file.tsv 2>&1 | tee $(basename $0 .sh).log\n\n"
    printf "Usage details:\n"
    printf "   1. run_mode: Either 'normal' for a single ATLAS assembly or 'coassembly' if the coassembly ATLAS extension was used (you'll know if you did this earlier; otherwise just use 'normal').\n"
    printf "   2. atlas_dir: Path to the base directory where ATLAS files were output.\n"
    printf "   3. assembly_sample_ID: Exact name of the assembly/coassembly that you desire to visualize bins from. If in coassembly mode, MUST be a sample within the 'coassembly' directory in the 'atlas_dir' folder.\n"
    printf "   4. output_dir: Path to save the anvio database and associated files to.\n"
    printf "   5. threads: number of threads to run.\n"
    printf "   6. mapping_guide_file: MUST be specified if you run in 'normal' mode but is not needed for 'coassembly' mode. A tab-separated file (with headers). Column 1 (sample_name) should be the ID of the metagenome you want to map onto the ATLAS assembly. Column 2 (raw_read_filepath) should be a comma-separated list of complete filepaths to the raw read files (e.g., R1, R2, se) for that metagenome - no spaces. You must supply either R1,R2 or R1,R2,se for the script to work.\n\n"
    printf "Additional usage notes:\n"
    printf "   * REQUIREMENTS: Must have anvio4 installed (ideally via conda). Must have the Anvi'o COG database set up via 'anvi-setup-ncbi-cogs'. Also, need R and an internet connection. See 'Installation' section above for the recommended method to get everything set up.\n"
    printf "   * If you use the 'coassembly' run_mode, this script is compataible with atlas-extensions version '1.0.22-coassembly-r3' output.\n"
    printf "\n* Example mapping_guide_file.tsv:\n"
    printf "sample_name\traw_read_filepaths\n"
    printf "NEIF_Jul17\t/home/data/metagenomes/NEIF/QC/NEIF_Jul17_QC_R1.fastq.gz,/home/data/metagenomes/NEIF/QC/NEIF_Jul17_QC_R2.fastq.gz,/home/data/metagenomes/NEIF/QC/NEIF_Jul17_QC_se.fastq.gz\n"
    printf "NEIF_Sep17\t/home/data/metagenomes/NEIF2/QC/NEIF_Sep17_QC_R1.fastq.gz,/home/data/metagenomes/NEIF2/QC/NEIF_Sep17_QC_R2.fastq.gz,/home/data/metagenomes/NEIF2/QC/NEIF_Sep17_QC_se.fastq.gz\n"
    printf "NW1_Jul17\t/home/data/metagenomes/NW1/QC/NW1_QC_R1.fastq.gz,/home/data/metagenomes/NW1/QC/NW1_QC_R2.fastq.gz,/home/data/metagenomes/NW1/QC/NW1_QC_se.fastq.gz\n\n"
    exit 1
fi

# Example
# atlas-to-anvi.sh "/Hippodrome/jmtsuji/180123_ELA111314_atlas_r1.0.22_plus_full" "CA-L227-2013" "/Hippodrome/jmtsuji/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/04_anvio/CA-L227-2013" 6 2>&1 | tee atlas-to-anvi.log

# Set variables from user input:
run_mode=$1
atlas_dir=$2
assembly_sample_ID=$3
output_dir=$4
threads=$5
if [ $# == 6 ]; then
	read_mapping_samples_table=$6
fi

function test_inputs {
	# TODO - add tests to check dependencies exist

	if [ ! -d ${atlas_dir} ]; then
	
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: atlas_dir directory '${atlas_dir}' does not exist. Is your atlas_dir correct? Exiting..."
		exit 1
	
	fi
	
	# Now check the run mode and then do sub-checks based on the run mode.
	if [ ${run_mode} = "coassembly" ]; then

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Running in coassembly mode."

		if [ ! -d ${atlas_dir}/coassembly ]; then
		
			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: directory '${atlas_dir}/coassembly' does not exist. Have you run atlas-coassembly? Exiting..."
			exit 1
		
		fi


		if [ ! -d ${atlas_dir}/coassembly/${assembly_sample_ID} ]; then
		
			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: directory '${atlas_dir}/coassembly/${assembly_sample_ID}' does not exist. Is your assembly_sample_ID '${assembly_sample_ID}' correct? Exiting..."
			exit 1
		
		fi

	elif [ ${run_mode} = "normal" ]; then

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Running in normal mode."

		if [ ! -d ${atlas_dir}/${assembly_sample_ID} ]; then
		
			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: directory '${atlas_dir}/${assembly_sample_ID}' does not exist. Is your assembly_sample_ID '${assembly_sample_ID}' correct? Exiting..."
			exit 1
		
		fi


		if [ ! -f ${read_mapping_samples_table} ]; then

			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: mapping_guide_file not found at '${read_mapping_samples_table}'. The guide file is required for normal mode. See help message (by running this script with no arguments) for the format of the guide file. Exiting..."
			exit 1

		fi

		# Roughly check the format of the guide file
		OFS=${IFS}
		IFS=" "
		expected_first_line="sample_name\traw_read_filepaths"
		actual_first_line=$(head -n 1 ${read_mapping_samples_table})
		if [ $(printf ${expected_first_line}) != ${actual_first_line} ]; then

			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: mapping_guide_file does not appear to be the correct format. The first line reads as '${actual_first_line}' instead of the expected '$(printf ${expected_first_line})'. The guide file is required for normal mode. See help message (by running this script with no arguments) for the format of the guide file. Exiting..."
			exit 1

		fi
		IFS=${OFS}

	else

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: run_mode must be either 'normal' or 'coassembly', but you entered '${run_mode}'. Exiting..."
		exit 1

	fi
	
	
	if [ -d ${output_dir} ]; then
	
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: directory '${output_dir}' already exists. Please delete before starting this script (if you don't care about the data in that folder), or choose a different output directory. Exiting..."
		exit 1
	
	fi

}

function assign_working_directory {
	# Determines the base directory where ATLAS files are stored based on whether the run mode is normal or coassembly

	if [ ${run_mode} = "normal" ]; then
		work_dir=${atlas_dir}/${assembly_sample_ID}

	elif [ ${run_mode} = "coassembly" ]; then
		work_dir=${atlas_dir}/coassembly/${assembly_sample_ID}

	else
		echo "Something's wrong. Sorry. Please contact code developers. Exiting..." # Really should never get to this point, because the run_mode was already checked in test_inputs above.
		exit 1
	fi

	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Assigned working directory to '${work_dir}'"

}

function clean_up_anvio_log {
	# Gets rid of the excessively long lines output by anvi'o -- although at risk of loss of information
	# input: filepath of anvio log as argument
	# output: will overwrite that log with a clean one

	# Grab input
	local anvio_logfile=$1

	# Simply log (get rid of lines with tons of spaces)
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Cleaning up anvio log '${anvio_logfile}'"
	grep -v "          " ${anvio_logfile} > ${anvio_logfile}_tmp
	mv ${anvio_logfile}_tmp ${anvio_logfile}

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
	python gff_parser.py ${work_dir}/annotation/prokka/${assembly_sample_ID}.gff \
					--gene-calls ${assembly_sample_ID}_gene_calls.txt --annotation ${assembly_sample_ID}_gene_annot.txt
	# Script gives no log info.

}

## TODO - add a different sanity check. Function in its current form does not work if genes without taxonomic assignment are at the end of the file...
function match_atlas_table_to_prokka_info {
	# TODO - move variables to local (within function) instead of global
	
	cd ${output_dir}/01b_import_atlas_table
	
	# Get the final gene entry number in the prokka file
	last_prokka_gene_ID=$(tail -n 1 ${output_dir}/01a_import_prokka/${assembly_sample_ID}_gene_calls.txt | cut -d $'\t' -f 1)
	
	# Get the final gene entry number of the altas taxonomy file
	last_atlas_gene_ID=$(tail -n 1 ${output_dir}/01b_import_atlas_table/${assembly_sample_ID}_gene_taxonomy.tsv | cut -d $'\t' -f 1)

	if [ ${last_prokka_gene_ID} -eq ${last_atlas_gene_ID} ]; then
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Prokka gene calls and ATLAS gene taxonomy annotations match in total numbers - good."
		
	elif [ ${last_prokka_gene_ID} -lt ${last_atlas_gene_ID} ]; then
		
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: WARNING: prokka gene calls and ATLAS gene taxonomy annotations do NOT match in total numbers. Final prokka gene call ID is '${last_prokka_gene_ID}', but final ATLAS taxonomy gene annotation ID is '${last_atlas_gene_ID}'. This has been observed before possibly as an error in prokka while generating the GFF file (final entries are cut out?)."

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: To correct error, deleting gene IDs in ATLAS gene taxonomy file from after ${last_prokka_gene_ID} until ${last_atlas_gene_ID}. Keeping old (complete) gene taxonomy file as '${assembly_sample_ID}_gene_taxonomy_ORIGINAL.tsv' but saving the new (truncated) gene taxonomy file over the original at '${assembly_sample_ID}_gene_taxonomy.tsv'."
		
		# Find the line number in the altas gene taxonomy file where the 'last_prokka_gene_ID' is
		atlas_tax_file="${output_dir}/01b_import_atlas_table/${assembly_sample_ID}_gene_taxonomy.tsv"
		
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
		mv ${assembly_sample_ID}_gene_taxonomy.tsv ${assembly_sample_ID}_gene_taxonomy_ORIGINAL.tsv
		head -n ${matching_atlas_tax_line} ${assembly_sample_ID}_gene_taxonomy_ORIGINAL.tsv > ${assembly_sample_ID}_gene_taxonomy.tsv
		
	elif [ ${last_prokka_gene_ID} -gt ${last_atlas_gene_ID} ]; then
	
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: prokka gene calls and ATLAS gene taxonomy annotations do NOT match in total numbers. Final prokka gene call ID is '${last_prokka_gene_ID}', but final ATLAS taxonomy gene annotation ID is '${last_atlas_gene_ID}'. The case where the ATLAS annotations are truncated compared to the prokka annotations has not been seen before. Something could be wrong with your ALTAS inputs, or this could be an unknown bug. Exiting..."
		exit 1
	
	fi
	
}

function export_atlas_info {

	cd ${output_dir}/01b_import_atlas_table

	# The annotations table has a different name after coassembly versus standard assembly
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Exporting information from the ATLAS annotations table (log: 01b_import_atlas_table/parse_atlas_table_for_anvio.log)"
	
	if [ ${run_mode} = "normal" ]; then
		annotations_filepath=${work_dir}/${assembly_sample_ID}_annotations.txt
	else
		annotations_filepath=${work_dir}/${assembly_sample_ID}_annotations_multi_mapping.txt
	fi

	if [ -f ${annotations_filepath} ]; then
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Found ATLAS annotations table at '${annotations_filepath}'."
	else
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Could not find ATLAS annotations table at '${annotations_filepath}'. Exiting..."
		exit 1
	fi

	# TODO -- Will replacing with the non multi-mapping version break anything? CHECK.
	parse_atlas_table_for_anvio.R -a ${annotations_filepath} \
					-t ${assembly_sample_ID}_gene_taxonomy.tsv -c ${assembly_sample_ID}_binning_results.tsv \
					-b ${assembly_sample_ID}_bins_info.tsv 2>&1 -d TRUE -@ ${threads} > parse_atlas_table_for_anvio.log

	# TODO - add a different sanity check. Function in its current form does not work if genes without taxonomic assignment are at the end of the file...
	## Deal with the differing length of the gene calls from the gff file versus the gene taxonomy
	# match_atlas_table_to_prokka_info

}

function generate_contig_database {

	#### 2. Generate contigs database
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Generating the contigs database (log: misc_logs/anvi-gen-contigs-database.log)"
	cd ${output_dir}
	anvi-gen-contigs-database -f ${work_dir}/${assembly_sample_ID}_contigs.fasta \
					-o ${assembly_sample_ID}_contigs.db -n ${assembly_sample_ID}_contigs_db \
					--external-gene-calls ${output_dir}/01a_import_prokka/${assembly_sample_ID}_gene_calls.txt \
					--ignore-internal-stop-codons --split-length -1 > misc_logs/anvi-gen-contigs-database.log 2>&1

	# Clean up log
	clean_up_anvio_log misc_logs/anvi-gen-contigs-database.log

	## Example table
	# gene_callers_id	contig	start	stop	direction	partial	source	version
	# 1	contig_01	1113	1677	f	0	program	v1.0
	## "The statement above means that the index of the first nucleotide in any contig should be 0. In other words, we start counting from 0, not from 1."

}

function add_hmm_annotations {
	
	#### Annotate with single-copy marker genes
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Annotating with single-copy marker genes (log: misc_logs/anvi-run-hmms.log)"
	cd ${output_dir}
	anvi-run-hmms -c ${assembly_sample_ID}_contigs.db --num-threads ${threads} \
					> misc_logs/anvi-run-hmms.log 2>&1
	# anvi-run-hmms -c ${assembly_sample_ID}_contigs.db --num-threads ${threads} --hmm-profile-dir [your_dir_for_custom_hmms]
	# anvi-display-contigs-stats ${assembly_sample_ID}_contigs.db # Will this work?
	
	# Clean up log
	clean_up_anvio_log misc_logs/anvi-run-hmms.log

}

function add_cog_annotations {
	
	# TODO - find a way to test whether or not setup is needed. Assumes already set up for now.
	
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Annotating with COGs (log: misc_logs/anvi-run-ncbi-cogs.log)"
	cd ${output_dir}
	mkdir tmp
	anvi-run-ncbi-cogs --num-threads ${threads} -c ${assembly_sample_ID}_contigs.db \
					--temporary-dir-path tmp > misc_logs/anvi-run-ncbi-cogs.log 2>&1 # --cog-data-dir ${cogs_data_dir}
	rm -rf tmp

	# Clean up log
	clean_up_anvio_log misc_logs/anvi-run-ncbi-cogs.log

}

function import_atlas_annotations {
	cd ${output_dir}

	# TODO - get past the error within anvi'o where this seems incompatible with COG annotations. Skipping for now.
	#### Import functional info
	#echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing functional annotations from prokka"
	#anvi-import-functions -c ${assembly_sample_ID}_contigs.db \
	#				-i ${output_dir}/01a_import_prokka/${assembly_sample_ID}_gene_annot.txt \
	#				2>&1 | tee misc_logs/anvi-import-functions.log

	## Example table
	# gene_callers_id	source	accession	function	e_value
	# 1	Pfam	PF01132	Elongation factor P (EF-P) OB domain	4e-23

	
	#### Import taxonomy info
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing taxonomic gene classifications from ATLAS (log: misc_logs/anvi-import-taxonomy.log)"
	anvi-import-taxonomy-for-genes -c ${assembly_sample_ID}_contigs.db \
                                        -i 01b_import_atlas_table/${assembly_sample_ID}_gene_taxonomy.tsv \
                                        -p default_matrix > misc_logs/anvi-import-taxonomy.log 2>&1

	## Example table
	# gene_callers_id	t_kingdom	t_phylum	t_class	t_order	t_family	t_genus	t_species
	# 1	 	 	 	 		 	Bacteroides fragilis
	# 2	 	 	 	 		 	Bacteroides fragilis

	# Clean up log
	clean_up_anvio_log misc_logs/anvi-import-taxonomy.log

}

function make_read_mapping_profiles_coassembly {

	# Get samples
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Mapping relative abundance profiles"
	cd ${output_dir}/02_multi_mapping
	sample_names=($(find ${work_dir}/multi_mapping -name "*.bam" -type f))

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
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: creating mapping profile (log: 02_multi_mapping/logs/anvi-profile_${sample_name_simple}.log)"
		anvi-profile -i ${sample_name_simple}.bam -c ${output_dir}/${assembly_sample_ID}_contigs.db \
						--output-dir ${sample_name_simple} --sample-name ${sample_name_simple} -T ${threads} \
						--min-contig-length 1000 > logs/anvi-profile_${sample_name_simple}.log 2>&1
		
		# Clean up log
		clean_up_anvio_log logs/anvi-profile_${sample_name_simple}.log

		# Remove indexed bam
		rm ${sample_name_simple}.bam ${sample_name_simple}.bam.bai
		
	done

}

function make_read_mapping_profiles_regular_assembly {

	# Get samples
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Mapping relative abundance profiles for samples provided by user"
	cd ${output_dir}/02_multi_mapping
	mkdir -p alignment_files

	# Load names of samples to map
	sample_names=($(cut -d $'\t' -f 1 ${read_mapping_samples_table} | tail -n +2))
	# Load a comma-separated list of the raw read inputs
	sample_filepaths=($(cut -d $'\t' -f 2 ${read_mapping_samples_table} | tail -n +2))
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${#sample_names[@]} samples provided."

	for i in $(seq 1 ${#sample_names[@]}); do
		j=$((${i}-1))
		sample_name=${sample_names[${j}]}
		sample_name_simple=$(sed s/'-'/'_'/g <<<${sample_name}) # Get rid of dashes
		sample_filepath=${sample_filepaths[${j}]}

		# Get the individual filepaths that are currently comma-separated. To do so, temporarily change the internal fields separator (IFS) to separate based on commas
		OFS=${IFS} # backup the standard internal fields separator
		IFS=","
		sample_filepaths_individual=(${sample_filepath})
		IFS=${OFS} # restore the old internal fields separator

		# Check files exist
		for path in ${sample_filepaths_individual[@]}; do
			if [ -f ${path} ]; then
				echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: found file '${path}'"
			else
				echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: ${sample_name}: file '${path}' not found. Exiting..."
				exit 1
			fi
		done

		# Assign generic variables for read mapping
		R1=${sample_filepaths_individual[0]}
		R2=${sample_filepaths_individual[1]}
		contigs="${work_dir}/${assembly_sample_ID}_contigs.fasta"
		outfile="alignment_files/${sample_name}_to_${assembly_sample_ID}.sam"
		logfile="alignment_files/${sample_name}_to_${assembly_sample_ID}_contig_coverage_stats.log"

		# Hard-code memory for now
		MEMORY=50

		# Do different read mapping based on whether two (R1, R2) or three file (R1, R2, se) paths were parsed out
		if [ ${#sample_filepaths_individual[@]} == 3 ]; then

			# Add single-ended reads
			se=${sample_filepaths_individual[2]}

			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: read mapping using ${#sample_filepaths_individual[@]} identified raw read files."
			bbwrap.sh nodisk=t ref=${contigs} in1=${R1},${se} in2=${R2},null perfectmode=t trimreaddescriptions=t outm=${outfile} threads=${threads} pairlen=1000 pairedonly=t mdtag=t xstag=fs nmtag=t sam=1.3 local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${MEMORY}G 2> ${logfile}

		elif [ ${#sample_filepaths_individual[@]} == 2 ]; then

			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: read mapping using ${#sample_filepaths_individual[@]} identified raw read files."
			bbwrap.sh nodisk=t ref=${contigs} in1=${R1} in2=${R2},null perfectmode=t trimreaddescriptions=t outm=${outfile} threads=${threads} pairlen=1000 pairedonly=t mdtag=t xstag=fs nmtag=t sam=1.3 local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${MEMORY}G 2> ${logfile}

		else
			
			# TODO - add this check to the start of the script so as to not disappoint the user.
			echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: ${sample_name}: user provided ${#sample_filepaths_individual[@]} raw read files for read mapping, but either 2 or 3 are needed. Exiting..."
			exit 1	
	
		fi

		# Convert SAM to BAM and index
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: Converting sam to bam"
		samtools view -@ ${threads} -u ${outfile} | samtools sort -@ ${threads} > ${outfile%.sam}.bam
		samtools index -b -@ ${threads} ${outfile%.sam}.bam
		rm ${outfile}
		
		# Generate profile
		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ${sample_name}: creating mapping profile (log: 02_multi_mapping/logs/anvi-profile_${sample_name_simple}.log)"
		anvi-profile -i ${outfile%.sam}.bam -c ${output_dir}/${assembly_sample_ID}_contigs.db \
						--output-dir ${sample_name_simple} --sample-name ${sample_name_simple} -T ${threads} \
						--min-contig-length 1000 > logs/anvi-profile_${sample_name_simple}.log 2>&1
		
		# Clean up log
		clean_up_anvio_log logs/anvi-profile_${sample_name_simple}.log

		# Remove indexed bam to save space
		rm ${outfile%.sam}.bam ${outfile%.sam}.bam.bai
		
	done

}

function merge_read_mapping_profiles {
	cd ${output_dir}

	# Check how many read mapping profiles are present
	profiles=($(find . -name "PROFILE.db"))
	num_profiles=${#profiles[@]}

	if [ ${num_profiles} = 1 ]; then

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Only one mapped sample; no need to merge mapped samples. Moving sample to '${assembly_sample_ID}_samples_merged'."
		profile=${profiles[0]}
		profile_folder=${profile%/*}
		mv ${profile_folder} ${assembly_sample_ID}_samples_merged

	elif [ ${num_profiles} -gt 1 ]; then

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Merging information from ${num_profiles} mapped samples into contig database (log: misc_logs/anvi-merge.log)."
		anvi-merge ${output_dir}/02_multi_mapping/*/PROFILE.db -o ${assembly_sample_ID}_samples_merged \
						-c ${assembly_sample_ID}_contigs.db --skip-concoct-binning -S metabat2 \
						> misc_logs/anvi-merge.log 2>&1
		# Consider '--enforce-hierarchical-clustering' to cluster even with > 25,000 contigs. But could take a long time...
		
		# Clean up log
		clean_up_anvio_log misc_logs/anvi-merge.log
		
	else

		echo "[$(date '+%y%m%d %H:%M:%S %Z')]: ERROR: found ${num_profiles} read mapping profiles to merge; this can't work for some reason. Exiting..."
		exit 1	

	fi

}

function import_custom_bins {

	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Importing bin info from ATLAS (log: misc_logs/anvi-import-collection.log)"
	cd ${output_dir}
	anvi-import-collection -p ${assembly_sample_ID}_samples_merged/PROFILE.db \
					-c ${assembly_sample_ID}_contigs.db -C "metabat2" --contigs-mode \
					--bins-info 01b_import_atlas_table/${assembly_sample_ID}_bins_info.tsv \
					01b_import_atlas_table/${assembly_sample_ID}_binning_results.tsv
					> misc_logs/anvi-import-collection.log 2>&1

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

	# Clean up log
	clean_up_anvio_log misc_logs/anvi-import-collection.log

}

function summarize {

	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: Summarizing anvio binning statistics (log: misc_logs/anvi-summarize.log)"
	cd ${output_dir}
	anvi-summarize -p ${assembly_sample_ID}_samples_merged/PROFILE.db \
					-c ${assembly_sample_ID}_contigs.db -C metabat2 \
					-o ${assembly_sample_ID}_summary --taxonomic-level t_genus \
					--init-gene-coverages > misc_logs/anvi-summarize.log 2>&1
	# Note: '--init-gene-coverages' takes a long time to run.
	# Note: consider '--quick-summary' for faster run with with more minimal output.

	# Clean up log
	clean_up_anvio_log misc_logs/anvi-summarize.log

}

function main {

	echo "Running $(basename $0), version ${script_version}, on $(date)."
	echo ""

	# Get date and time of start
	start_time=$(date)
	
	# Report settings
	echo "Run settings:"
	echo "atlas_dir: ${atlas_dir}"
	echo "assembly_sample_ID: ${assembly_sample_ID}"
	echo "output_dir: ${output_dir}"
	echo "threads: ${threads}"
	echo ""

	test_inputs
	assign_working_directory
	
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
	if [ ${run_mode} == "coassembly" ]; then
		make_read_mapping_profiles_coassembly
	elif [ ${run_mode} == "normal" ]; then
		make_read_mapping_profiles_regular_assembly
	fi

	merge_read_mapping_profiles

	# Import bins
	import_custom_bins
	
	# Generate summary
	summarize
	
	end_time=$(date)

	# Exit and report anvi-interative instructions to user
	echo "[$(date '+%y%m%d %H:%M:%S %Z')]: $(basename $0): complete. Started at ${start_time} and finished at ${end_time}.\n\n"
	printf "# To visualize, please run:\n"
	printf "cd ${output_dir}\n"
	printf "# For overview of all bins:\n"
	printf "anvi-interactive -p ${assembly_sample_ID}_samples_merged/PROFILE.db -c ${assembly_sample_ID}_contigs.db -C metabat2 --server-only -P 8080\n"
	printf "# For refining a single bin:\n"
	printf "bin_name=[name_of_interest]\n"
	printf "anvi-refine -p ${assembly_sample_ID}_samples_merged/PROFILE.db -c ${assembly_sample_ID}_contigs.db -C metabat2 -b \${bin_name} --taxonomic-level t_genus --title \${bin_name} -P 8080 --server-only\n"
	printf "# You can vary the taxonomic level visualized here -- e.g., t_genus, t_family, and so on."
	printf "# To see the anvi'o plot generated by the above visualization commands, go onto your web browser after running and go to [your_server_url]:8080 to view. E.g., neufeldserver.uwaterloo.ca:8080. I'm assuming that you want to stream the results from a remote server here (rather than just running anvi'o on your own laptop).\n\n"
	printf "# Then, to export the refined bin info, run:\n"
	printf "anvi-summarize -p ${assembly_sample_ID}_samples_merged/PROFILE.db -c ${assembly_sample_ID}_contigs.db -C metabat2 -o ${assembly_sample_ID}_summary_refined --taxonomic-level t_family --init-gene-coverages 2>&1 | tee misc_logs/anvi-summarize-refined.log\n\n"

}

main

