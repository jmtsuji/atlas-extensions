#!/bin/bash
# co-assembly.sh
# Copyright Jackson M. Tsuji, 2018
# Neufeld lab, University of Waterloo, Canada
# Created Jan. 29th, 2018
# Description: Re-arranges the ATLAS pipeline (in crude form, for now) to perform co-assembly on specific samples.
# REQUIRES that ATLAS (e.g., version 1.0.22) has been run normally on all samples first

# Basic script stuff (from Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 397):
set -e
set -u
set -o pipefail

script_version="1.0.22r1" # to match ATLAS version it is designed to work with

# If no input is provided, exit out and provide help
if [ $# == 0 ]
    then
    printf "$(basename $0): ______.\n"
    printf "Version: ${script_version}\n"
    printf "Contact Jackson M. Tsuji (jackson.tsuji@uwaterloo.ca) for bug reports or feature requests.\n\n"
    printf "Usage: $(basename $0) output_dir database_dir config_filepath coassembly_guide_filepath threads\n\n"
    printf "Usage details:\n"
    printf "1. output_dir: Path to the base directory where ATLAS files were output.\n"
    printf "2. database_dir: Path to the ATLAS database directory.\n"
    printf "3. config_filepath: Text list of the ATLAS sample names of samples that should be combined for co-assembly.\n"
    printf "4. coassembly_guide_filepath: TSV file with three columns: coassembly_name (names of coassembly runs); assembly_samples (comma-separated names of previously run ATLAS samples to coassemble); read_mapping_samples (comma-separated names of previously run ATLAS samples to read map for binning).\n"
    printf "5. threads: number of threads to run.\n\n"
    exit 1
fi
# Using printf: http://stackoverflow.com/a/8467449 (accessed Feb 21, 2017)
# Test for empty variable: Bioinformatics Data Skills Ch. 12 pg 403-404, and http://www.tldp.org/LDP/Bash-Beginners-Guide/html/sect_09_07.html and http://stackoverflow.com/a/2428006 (both accessed Feb 21, 2017)

# Set variables from user input:
OUTPUT_DIR=$1
DATABASE_DIR=$2
CONFIG_FILEPATH=$3
COASSEMBLY_GUIDE_FILEPATH=$4
THREADS=14 # TODO Eventually need to set this to user input

# Constants
ASSEMBLY_BINARIES="${OUTPUT_DIR}/.snakemake/conda/ccef13c9/bin"
BINNING_BINARIES="${OUTPUT_DIR}/.snakemake/conda/a6a01fdc/bin"


function test_inputs {
	# Description: tests that provided folders and files exist in the proper configuration
	# Params: No input needed; relies on globally assigned variables above
	
	if [ ! -d $OUTPUT_DIR ]; then
		echo "ERROR: Cannot find output directory at '${OUTPUT_DIR}'. Exiting..."
		exit 1
	fi
	
	if [ ! -d $DATABASE_DIR ]; then
		echo "ERROR: Cannot find database directory at '${DATABASE_DIR}'. Exiting..."
		exit 1
	fi
	
	if [ ! -f $COASSEMBLY_GUIDE_FILEPATH ]; then
		echo "ERROR: Cannot find coassembly guide file at ${COASSEMBLY_GUIDE_FILEPATH}. Exiting..."
		exit 1
	fi
	
	if [ ! -d $ASSEMBLY_BINARIES ]; then
		echo "ERROR: Cannot find conda binaries for assembly in ${ASSEMBLY_BINARIES}. Have you run ATLAS normally yet? Exiting..."
		exit 1
	fi

	if [ ! -d $BINNING_BINARIES ]; then
		echo "ERROR: Cannot find conda binaries for genome binning in ${BINNING_BINARIES}. Have you run ATLAS normally yet? Exiting..."
		exit 1
	fi
		
}

function get_coassembly_params {
	# Description: parses the coassembly guide file to get coassembly names, files to assembly, files to map. Checks they exist.
	# GLOBAL params: COASSEMBLY_GUIDE_FILEPATH, OUTPUT_DIRECTORY
	# Local params: none
	# Return: coassembly_names (array); assembly_samples (array); read_mapping_samples (array)
	
	# Check headers of COASSEMBLY_GUIDE
	local coassembly_headers=($(head -n 1 ${COASSEMBLY_GUIDE_FILEPATH}))
	
	if [ ${#coassembly_headers[@]} != 3 ]; then
		# TODO improve this error message.
		echo "ERROR: coassembly guide file has incorrect number of columns. Looking for 3 but found ${#coassembly_headers[@]}. Job terminating."
		exit 1
	elif [ ${coassembly_headers[0]} != "coassembly_name" -o  ${coassembly_headers[1]} != "assembly_samples" \
		-o ${coassembly_headers[2]} != "read_mapping_samples" ]; then
		# TODO improve this error message.
		echo "ERROR: coassembly guide file has incorrect headers. Job terminating."
		exit 1
	fi
	
	# Get COASSEMBLY_GUIDE params
	coassembly_names=($(tail -n +2 ${COASSEMBLY_GUIDE_FILEPATH} | cut -d $'\t' -f 1))
	assembly_samples=($(tail -n +2 ${COASSEMBLY_GUIDE_FILEPATH} | cut -d $'\t' -f 2))
	read_mapping_samples=($(tail -n +2 ${COASSEMBLY_GUIDE_FILEPATH} | cut -d $'\t' -f 3))
	
	# Test COASSEMBLY_GUIDE params are okay
	check_coassembly_params
	
}

function check_coassembly_params {
	# Description: checks that each specified sample has been run through ATLAS and that coassembly directories do not already exist (makes them at the same time)
	# GLOBAL params: OUTPUT_DIRECTORY; coassembly_names (array); assembly_samples (array); read_mapping_samples (array)
	# Local params: none
	# Return: none

	check_coassembly_dirs
	check_samples assembly
	check_samples mapping

}

function check_samples {
	# Description: tests if each input assembly/mapping sample exists
	# GLOBAL params: OUTPUT_DIR, assembly_samples OR read_mapping_samples (arrays containing strings of comma-separated values)
	# Local params: sample_type (string; either 'assembly' or 'mapping')
	# Return: none

	local sample_type=$1

	if [ $sample_type == "assembly" ]; then
		local samples=(${assembly_samples[@]})
	elif [ $sample_type == "mapping" ]; then
		local samples=(${read_mapping_samples[@]})
	else
		echo "ERROR: either 'assembly' or 'mapping' should be passed to check_samples. Instead, got '${sample_type}'. Job terminating."
		exit 1
	fi

	# Each coassembly_name has several associated assembly_samples separated by commas. Examine one coassembly_name at a time.
	for coassembly_sample in ${samples[@]}; do

		# Temporarily change the internal fields separator (IFS) to parse comma separators. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
		local OFS="$IFS"
		IFS=,
		
		# Get names of individual samples provided for that coassembly_name
		local assembly_sample_IDs=(${coassembly_sample})
		
		# Iteratively check each assembly_sample to see if it exists and exit if not.
		for sample_ID in ${assembly_sample_IDs[@]}; do
			local atlas_path="${OUTPUT_DIR}/${sample_ID}"
			local qc_filepath="${atlas_path}/sequence_quality_control"
			
			if [ ! -d ${atlas_path} ]; then
				echo "ERROR: ${sample_ID}: Cannot find path for assembly sample. Looking in '${atlas_path}'. Exiting..."
				exit 1
			elif [ ! -f ${qc_filepath}/${sample_ID}_QC_se.fastq.gz ]; then
				# TODO: confirm that se reads will be present both for paired end and R1/R2. This test might not be sufficient right now.
				echo "ERROR: ${sample_ID}: Cannot find finished QC files '${sample_ID}_QC_se.fastq.gz' in '${qc_filepath}'. Exiting..."
				exit 1 
			fi
			
		done
	
		# Fix the IFS
		IFS="$OFS"
	
	done
		
}

function check_coassembly_dirs {
	# Description: tests if coassembly directories exist (they should not!) and makes directories otherwise.
	# GLOBAL params: OUTPUT_DIR, coassembly_names (array)
	# Local params: none
	# Return: none
	
	mkdir -p ${OUTPUT_DIR}/coassembly
	
	for name in ${coassembly_names[@]}; do
		local coassembly_dir="${OUTPUT_DIR}/coassembly/${name}"
		
		if [ -d ${coassembly_dir} ]; then
			echo "ERROR: found pre-existing output co-assembly directory ${coassembly_dir}. Please delete before starting. Job terminating."
			exit 1
		fi

	done

}


function check_yaml {
	# Description: If the provided .yaml filepath does not yet exist, generates yaml file for user based on coassembly samples, then exits for user to modify. If the .yaml filepath does exist, then moves forward with the analysis.
	# GLOBAL params: OUTPUT_DIR; DATABASE_DIR; CONFIG_FILEPATH
	# Local params: none
	# Return: writes .yaml to disk if not already existing

	if [ ! -f ${CONFIG_FILEPATH} ]; then
		
		build_yaml
		
		echo "Created .yaml configuration file for the coassembly run at '${CONFIG_FILEPATH}'. \
			Please edit this using a text editor to tweak any desired settings. \
			Once ready to go, specify the edited version as the config_coassembly.yaml argument \
			into this script (keeping all other arguments the same as the current run). Exiting for now..."
		# TODO add example code for what the new run should look like, to guide the user.
		
		exit 1

	elif [ -f ${CONFIG_FILEPATH} ]; then

		# Perform sanity check on provided .yaml to make sure it's actually for the coassemblies and not a mistake (i.e. original ATLAS file)
		if ! grep -q ${coassembly_names[0]} ${CONFIG_FILEPATH}; then
			echo "ERROR: provided configuration file '${CONFIG_FILEPATH}' does not contain coassembly sample '${coassembly_names[0]}'. Are you sure you provided the right configuration file? Exiting..."
			exit 1
		fi
		
		echo "Configuration (.yaml) file for coassemblies provided. Not generating a new one."
	
	else

		echo "Something is wrong with the provided path for the .yaml configuration file. Job terminating."
		exit 1

	fi

}


function build_yaml {
	# Description: Generates yaml file for user based on coassembly samples.
	# GLOBAL params: OUTPUT_DIR; DATABASE_DIR; CONFIG_FILEPATH
	# Local params: none
	# Return: writes .yaml to disk
	
	local temp_coassembly_dir="${OUTPUT_DIR}/coassembly/inputs_tmp"
	
	# Check if temp directory already exists (so that the script does not erase it later):
	if [ ! -d ${temp_coassembly_dir} ]; then
		local erase_temp_coassembly_dir="TRUE"
	fi
	# TODO perform this check later??
	
	# Make temporary directory for fake samples from which to generate the config file
	mkdir -p ${temp_coassembly_dir}
	
	# Make fake samples
	for name in ${coassembly_names[@]}; do
		touch ${temp_coassembly_dir}/${name}_R1.fastq.gz
		touch ${temp_coassembly_dir}/${name}_R2.fastq.gz
	done
	
	# Generate config file from fake samples
	atlas make-config --database-dir ${DATABASE_DIR} ${CONFIG_FILEPATH} ${temp_coassembly_dir}
	
}


function concatenate_pre_assembly {
	# Description: concatenates decontaminated read files from the individual ATLAS runs (according to assembly_names.list) in preparation for co-assembly
	# GLOBAL Params: 'assembly_names' (array) from 'get_assembly_samples'
	# Local params: none
	# Return: writes files to disk (.fastq.gz)

	echo "Concatenating pre-assembly files..."

	# For each coassembly that will be performed:
	for i in $(seq 1 ${#coassembly_names[@]}); do
		# Set counter to be based on zero, not one
		local j=$((${i}-1))

		local coassembly_name=${coassembly_names[${j}]}

		# Make output directory for that coassembly sample
		local coassembly_dir="${OUTPUT_DIR}/coassembly/${coassembly_name}"
		mkdir -p ${coassembly_dir}/sequence_quality_control

		# Temporarily change the internal fields separator (IFS) to parse comma separators. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
		local OFS="$IFS"
		IFS=,

		# Get names of individual samples provided for that coassembly_name
		local assembly_sample_IDs=(${assembly_samples[${j}]})

		echo "${coassembly_name}: ${assembly_sample_IDs[@]}"

		# Concatenate the reads from each individual sample
		for sample in ${assembly_sample_IDs[@]}; do

			cat ${OUTPUT_DIR}/${sample}/sequence_quality_control/${sample}_QC_R1.fastq.gz >> ${coassembly_dir}/sequence_quality_control/${coassembly_name}_QC_R1.fastq.gz
			cat ${OUTPUT_DIR}/${sample}/sequence_quality_control/${sample}_QC_R2.fastq.gz >> ${coassembly_dir}/sequence_quality_control/${coassembly_name}_QC_R2.fastq.gz
			cat ${OUTPUT_DIR}/${sample}/sequence_quality_control/${sample}_QC_se.fastq.gz >> ${coassembly_dir}/sequence_quality_control/${coassembly_name}_QC_se.fastq.gz		

		done

	# Fix the IFS
	IFS="$OFS"

	echo ""
	done
}

function make_atlas_dirs {
	# Description: forces ATLAS to think that run QC has already been performed on the provided samples
	
	local coassembly_dir="${OUTPUT_DIR}/coassembly"
	
	# Make fake samples
	for sample in ${coassembly_names[@]}; do
		mkdir -p ${coassembly_dir}/${sample}/logs ${coassembly_dir}/ref/genome/1 ${coassembly_dir}/logs
		
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_raw_R1.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_raw_R2.fastq.gz
		
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_deduplicated_R1.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_deduplicated_R2.fastq.gz
		
		
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_filtered_R1.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_filtered_R2.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_filtered_se.fastq.gz
		touch ${coassembly_dir}/${sample}/logs/${sample}_quality_filtering_stats.txt
		
		touch ${coassembly_dir}/ref/genome/1/summary.txt
		
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_clean_R1.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_clean_R2.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_clean_se.fastq.gz
		
		mkdir -p ${coassembly_dir}/${sample}/sequence_quality_control/contaminants
		touch ${coassembly_dir}/${sample}/sequence_quality_control/contaminants/PhiX_R1.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/contaminants/PhiX_R2.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/contaminants/PhiX_se.fastq.gz
		
		touch ${coassembly_dir}/${sample}/sequence_quality_control/contaminants/rRNA_R1.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/contaminants/rRNA_R2.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/contaminants/rRNA_se.fastq.gz
		
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_decontamination_reference_stats.txt
		
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_QC_R1.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_QC_R2.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_QC_se.fastq.gz
		
	done

}


# TODO: UNTESTED
function run_atlas {
	# Description: runs standard ATLAS pipeline, starting from assembly, for coassembly files
	# TODO: make sure this strategy actually works!!!
	
	local coassembly_dir="${OUTPUT_DIR}/coassembly"
	local log_code=$(date '+%y%m%d_%H%M')
	
	# Define the step to force ATLAS to start from and end at
	local start_step="normalize_coverage_across_kmers"
	local end_step="sort_munged_blast_hits"
	
	# See what all steps of ATLAS would be without actually running ATLAS
	atlas assemble --jobs ${THREADS} --out-dir ${coassembly_dir} ${CONFIG_FILEPATH} --force ${start_step} --until ${end_step} --dryrun > ${coassembly_dir}/atlas_run_steps_${log_code}.log
	
	# Run ATLAS
	atlas assemble --jobs ${THREADS} --out-dir ${coassembly_dir} ${CONFIG_FILEPATH} --force ${start_step} --until ${end_step} 2>&1 | tee ${coassembly_dir}/atlas_run_${log_code}.log

}


# TODO: Unfinished
function main {
	echo "Running $(basename $0), version ${script_version}, on $(date)."
	echo ""
	
	test_inputs
	
	# Get inputs
	get_coassembly_params

	# Generate yaml and exit early, or continue on if yaml exists
	check_yaml

	# Get date and time of start
	start_time=$(date)

	concatenate_pre_assembly
	make_atlas_dirs
	run_atlas

	# TODO functions not even started yet...
	read_map_to_coassemblies
	bin_coassemblies
	# TODO can I run ATLAS on these bins for annotation in place of the standard MaxBin ones? Move the MaxBin ones elsewhere? (Or don't even generate them -- just have ATLAS stop before this point, then resume it)

	finish_atlas # to finish annotation on bins as stated above.
		

	end_time=$(date)	
	echo ""
	echo "Coassembly finished. Output saved in ${OUTPUT_DIR}."
	echo "Started at ${start_time} and finished at ${end_time}."
	echo ""

	echo "$(basename $0): finished."
	echo ""

}

main

