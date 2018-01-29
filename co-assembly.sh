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
    printf "\n$(basename $0): finds and replaces selected text in a file.\n"
    printf "Version: ${script_version}\n"
    printf "Contact Jackson M. Tsuji (jackson.tsuji@uwaterloo.ca) for bug reports or feature requests.\n\n"
    printf "Usage: $(basename $0) output_dir database_dir assembly_samples.list mapping_samples.list\n\n"
    printf "Usage details:\n"
    printf "- output_dir: Path to the base directory where ATLAS files were output.\n"
    printf "- database_dir: Path to the ATLAS database directory.\n"
    printf "- assembly_samples.list: Text list of the ATLAS sample names of samples that should be combined for co-assembly.\n"
    printf "- mapping_samples.list: Text list of the ATLAS sample names of samples that should be read mapped onto the co-assembly before binning.\n\n"
    exit 1
fi
# Using printf: http://stackoverflow.com/a/8467449 (accessed Feb 21, 2017)
# Test for empty variable: Bioinformatics Data Skills Ch. 12 pg 403-404, and http://www.tldp.org/LDP/Bash-Beginners-Guide/html/sect_09_07.html and http://stackoverflow.com/a/2428006 (both accessed Feb 21, 2017)

# Get date and time
DATETIME=$(date '+%y%m%d_%H%M')

# Set variables from user input:
OUTPUT_DIR=$1
DATABASE_DIR=$2
ASSEMBLY_SAMPLES=$3
MAPPING_SAMPLES=$4

# Constants
ASSEMBLY_BINARIES="${OUTPUT_DIR}/.snakemake/conda/ccef13c9/bin"
BINNING_BINARIES="${OUTPUT_DIR}/.snakemake/conda/a6a01fdc/bin"
CO_ASSEMBLY_NAME="co_assembly_${DATETIME}"
CO_ASSEMBLY_DIR="${OUTPUT_DIR}/${CO_ASSEMBLY_NAME}"

function test_inputs {
	# Description: tests that provided folders and files exist in the proper configuration
	# Params: No input needed; relies on globally assigned variables above
	
	if [ ! -d $OUTPUT_DIR ]; then
		echo "ERROR: Cannot find provided output directory ${OUTPUT_DIR}. Exiting..."
		quit 1
	fi
	
	if [ ! -d $DATABASE_DIR ]; then
		echo "ERROR: Cannot find provided database directory ${DATABASE_DIR}. Exiting..."
		quit 1
	fi
	
	if [ ! -f $ASSEMBLY_SAMPLES ]; then
		echo "ERROR: Cannot find provided assembly sample names list ${ASSEMBLY_SAMPLES}. Exiting..."
		quit 1
	fi
	
	if [ ! -f $MAPPING_SAMPLES ]; then
		echo "ERROR: Cannot find provided mapping sample names list ${MAPPING_SAMPLES}. Exiting..."
		quit 1
	fi
	
	if [ ! -d $ASSEMBLY_BINARIES ]; then
		echo "ERROR: Cannot find conda binaries for assembly in ${ASSEMBLY_BINARIES}. Have you run ATLAS normally yet? Exiting..."
		quit 1
	fi

	if [ ! -d $BINNING_BINARIES ]; then
		echo "ERROR: Cannot find conda binaries for genome binning in ${BINNING_BINARIES}. Have you run ATLAS normally yet? Exiting..."
		quit 1
	fi
	
	# Make co-assembly directory
	if [ ! -d ${CO_ASSEMBLY_DIR} ]; then
		mkdir ${CO_ASSEMBLY_DIR}
		echo "Results will be stored in ${CO_ASSEMBLY_DIR}"
		echo ""
	else
		echo "ERROR: found pre-existing output co-assembly directory ${CO_ASSEMBLY_DIR}. Please delete before starting. Exiting..."
		quit 1
	fi
		
}

function get_assembly_samples {
	# Description: parses the input assembly sample names list into an array and checks they exist
	# Params: list of assembly sample names (.list)
	# Return: assembly_names (array of assembly_names)
	
	local assembly_list=$1
	
	assembly_names=($(cat ${assembly_list}))
	
	for name in ${assembly_names[@]}; do
		atlas_path="${OUTPUT_DIR}/${name}"
		
		if [ ! -d ${atlas_path} ]; then
			echo "ERROR: Cannot find path for required assembly sample ${atlas_path}. Exiting..."
			quit 1
		fi
	done
	
}

function get_mapping_samples {
	# Description: parses the input mapping sample names list into an array and checks they exist
	# Params: list of mapping sample names (.list)
	# Return: mapping_names (array of assembly_names)
	
	local mapping_list=$1
	
	mapping_names=($(cat ${mapping_list}))
	
	for name in ${mapping_names[@]}; do
		atlas_path="${OUTPUT_DIR}/${name}"
		
		if [ ! -d ${atlas_path} ]; then
			echo "ERROR: Cannot find path for required mapping sample ${atlas_path}. Exiting..."
			quit 1
		fi
	done
	
}

function concatenate_pre_assembly {
	# Description: concatenates decontaminated read files from the individual ATLAS runs (according to assembly_names.list) in preparation for co-assembly
	# Params: requires that 'assembly_names' array is globally available from 'get_assembly_samples'

	echo "Concatenating pre-assembly files..."

	local out_dir=${CO_ASSEMBLY_DIR}/sequence_quality_control
	mkdir -p ${out_dir}

	for name in ${assembly_names[@]}; do
		cat ${OUTPUT_DIR}/${name}/sequence_quality_control/${name}_QC_R1.fastq.gz >> ${out_dir}/${CO_ASSEMBLY_NAME}_QC_R1.fastq.gz
		cat ${OUTPUT_DIR}/${name}/sequence_quality_control/${name}_QC_R2.fastq.gz >> ${out_dir}/${CO_ASSEMBLY_NAME}_QC_R2.fastq.gz
		cat ${OUTPUT_DIR}/${name}/sequence_quality_control/${name}_QC_se.fastq.gz >> ${out_dir}/${CO_ASSEMBLY_NAME}_QC_se.fastq.gz		
	done

}

function normalize_coverage_across_kmers {
	
	mkdir -p ${CO_ASSEMBLY_DIR}/assembly/reads
	mkdir -p ${CO_ASSEMBLY_DIR}/logs
	
	if [ in=${CO_ASSEMBLY_DIR}/sequence_quality_control/${CO_ASSEMBLY_NAME}_QC_se.fastq.gz != "null" ]; then
		bbnorm.sh in=${CO_ASSEMBLY_DIR}/sequence_quality_control/${CO_ASSEMBLY_NAME}_QC_se.fastq.gz extra=${CO_ASSEMBLY_DIR}/sequence_quality_control/${CO_ASSEMBLY_NAME}_QC_R1.fastq.gz,${CO_ASSEMBLY_DIR}/sequence_quality_control/${CO_ASSEMBLY_NAME}_QC_R2.fastq.gz out=${CO_ASSEMBLY_DIR}/assembly/reads/normalized_se.fastq.gz k=21 t=100 interleaved=f minkmers=15 prefilter=t threads=14 -Xmx65G 2> ${CO_ASSEMBLY_DIR}/logs/${CO_ASSEMBLY_NAME}_normalization.log
	fi



# NOT DONE...
	if [ t = "t" ]; then
		bbnorm.sh in=NE1-sub100k/sequence_quality_control/NE1-sub100k_QC_R1.fastq.gz in2=NE1-sub100k/sequence_quality_control/NE1-sub100k_QC_R2.fastq.gz extra=NE1-sub100k/sequence_quality_control/NE1-sub100k_QC_se.fastq.gz out=NE1-sub100k/assembly/reads/normalized_R1.fastq.gz out2=NE1-sub100k/assembly/reads/normalized_R2.fastq.gz k=21 t=100 interleaved=f minkmers=15 prefilter=t threads=14 -Xmx32G 2>> NE1-sub100k/logs/NE1-sub100k_normalization.log
	fi
}

function main {
	echo "Running $(basename $0), version ${script_version}, on $(date)."
	echo ""
	
	# Get inputs
	test_inputs
	get_assembly_samples # Now in variable 'assembly_names'
	get_mapping_samples # Now in variable 'mapping_names'
	echo "Running co-assembly for ${#assembly_names[@]} samples and mapping with ${#mapping_names[@]} samples."
	echo ""
	
	concatenate_pre_assembly

}


# Test that the ATLAS conda binaries are where they should be


# By default, set working directory to present working directory (pwd)
work_dir=$(pwd)

echo "Running $(basename $0), version $script_version."
echo ""
echo "Replacing items in file $(basename ${input_name}), based on the list provided in $(basename ${replacement_info})."


cd "${work_dir}"

# Temporarily change the internal fields separator (IFS) so that whitespaces in find/replace scheme do not create new entries. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
OFS="$IFS"
IFS=$'\n'

# Get lists of old and new filenames and store as arrays
names_old=($(tail -n +2 "${replacement_info}" | cut -d $'\t' -f 1))
names_new=($(tail -n +2 "${replacement_info}" | cut -d $'\t' -f 2))
# Got help from http://stackoverflow.com/questions/2961673/find-missing-argument-to-exec, post by Marian on June 2, 2010; accessed May 13, 2016

# Fix the IFS
IFS="$OFS"

# Get number of names
num_names=${#names_old[@]}

# Confirm number of names is equal between lists
if [ ${#names_old[@]} != ${#names_new[@]} ]
then
    echo "Error: ${#names_old[@]} entries in old names column and ${#names_new[@]} in new names column. These do not match. This could possibly be due to special characers in one of the columns, which this has limited support for. Exiting..."
    exit 1
else echo "Identified ${num_names} items for replacement."
fi

echo ""

# Make a copy of the tree file for multiple in-place edits
cp "${input_name}" "${output_name}"

for i in $(seq 1 ${num_names})
do
# Get names
    old_name="${names_old[i-1]}"
    new_name="${names_new[i-1]}"
    echo "${old_name} --> ${new_name}"
    sed -i -e "s/${old_name}/${new_name}/g" $output_name
done
# "for" loop idea from Bioinformatics Data Skills (Vince Buffalo), Ch. 12
# sed help from http://unix.stackexchange.com/a/159369 (main) and http://stackoverflow.com/a/15236526 (debugging), accessed Feb. 22, 2017

# Remove backup file created by sed once done
rm ${output_name}-e


echo ""
echo "Replacing finished. Output saved as $(basename ${output_name})."
echo ""

echo "$(basename $0): finished."
echo ""
