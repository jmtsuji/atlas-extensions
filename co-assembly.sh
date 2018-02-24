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
		
		printf "\nCreated .yaml configuration file for the coassembly run at '${CONFIG_FILEPATH}'.\n"
		printf "Please edit this using a text editor to tweak any desired settings.\n"
		printf "Once ready to go, specify the edited version as the config_coassembly.yaml argument"
		printf "into this script (keeping all other arguments the same as the current run). Exiting for now...\n\n"
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
	# GLOBAL Params: 'assembly_names' (array) from 'get_coassembly_params'
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
	# GLOBAL Params: OUTPUT_DIR; 'coassembly_names' (array) from 'get_coassembly_params'
	# Local params: none
	# Return: writes empty files/directories to disk
	
	local coassembly_dir="${OUTPUT_DIR}/coassembly"
	
	# Make fake samples, round 1
	for sample in ${coassembly_names[@]}; do
		mkdir -p ${coassembly_dir}/${sample}/logs ${coassembly_dir}/${sample}/sequence_quality_control/read_stats
		mkdir -p ${coassembly_dir}/ref/genome/1 ${coassembly_dir}/ref/index/1
		mkdir -p ${coassembly_dir}/logs/benchmarks
		
		touch ${coassembly_dir}/${sample}/logs/${sample}_init.log
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_raw_R1.fastq.gz
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_raw_R2.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/raw_read_counts.tsv
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/raw.zip
		
		mkdir -p ${coassembly_dir}/logs/benchmarks/deduplicate
		touch ${coassembly_dir}/logs/benchmarks/deduplicate/${sample}.txt
		touch ${coassembly_dir}/${sample}/logs/${sample}_deduplicate.log
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_deduplicated_R1.fastq.gz
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_deduplicated_R2.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/deduplicated_read_counts.tsv
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/deduplicated.zip
		
		mkdir -p ${coassembly_dir}/logs/benchmarks/quality_filter
		touch ${coassembly_dir}/logs/benchmarks/quality_filter/${sample}.txt
		touch ${coassembly_dir}/${sample}/logs/${sample}_quality_filter.log
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_filtered_R1.fastq.gz
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_filtered_R2.fastq.gz
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_filtered_se.fastq.gz
		touch ${coassembly_dir}/${sample}/logs/${sample}_quality_filtering_stats.txt
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/filtered_read_counts.tsv
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/filtered.zip
	
	done

	# Perform common step
	touch ${coassembly_dir}/ref/index/1/chr1_index_k13_c4_b1.block
	touch ${coassembly_dir}/ref/index/1/chr1_index_k13_c4_b1.block2.gz
	touch ${coassembly_dir}/ref/genome/1/summary.txt
	touch ${coassembly_dir}/ref/genome/1/info.txt
	touch ${coassembly_dir}/ref/genome/1/namelist.txt
	touch ${coassembly_dir}/ref/genome/1/merged_ref_5214708240339.fa.gz
	touch ${coassembly_dir}/ref/genome/1/reflist.txt
	touch ${coassembly_dir}/ref/genome/1/scaffolds.txt.gz
	touch ${coassembly_dir}/ref/genome/1/chr1.chrom.gz
	touch ${coassembly_dir}/logs/build_decontamination_db.log
	
	# Make fake samples, round 2
	for sample in ${coassembly_names[@]}; do
		mkdir -p ${coassembly_dir}/logs/benchmarks/decontamination
		touch ${coassembly_dir}/logs/benchmarks/decontamination/${sample}.txt
		touch ${coassembly_dir}/${sample}/logs/${sample}_decontamination.log
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_clean_R1.fastq.gz
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_clean_R2.fastq.gz
		#touch ${coassembly_dir}/${sample}/sequence_quality_control/${sample}_clean_se.fastq.gz
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/clean_read_counts.tsv
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/clean.zip
		
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
		touch ${coassembly_dir}/${sample}/logs/read_stats.log
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/QC_read_counts.tsv
		touch ${coassembly_dir}/${sample}/sequence_quality_control/read_stats/QC.zip
		
	done

}


function run_atlas_assemble {
	# Description: runs standard ATLAS pipeline, starting from assembly and finishing before read mapping, for coassembly files
	# GLOBAL Params: OUTPUT_DIR; THREADS; CONFIG_FILEPATH
	# Local params: none
	# Return: writes files/directories to disk

	local coassembly_dir="${OUTPUT_DIR}/coassembly"
	local log_code=$(date '+%y%m%d_%H%M')
	
	# Define the step to force ATLAS to start from and end at
	# local end_step="sort_munged_blast_hits"
	local end_step="make_maxbin_abundance_file"
	
	# See what all steps of ATLAS would be without actually running ATLAS
	atlas assemble --jobs ${THREADS} --out-dir ${coassembly_dir} ${CONFIG_FILEPATH} --until ${end_step} --dryrun > ${coassembly_dir}/atlas_run_steps_${log_code}.log 2>&1
	
	# Run ATLAS
	atlas assemble --jobs ${THREADS} --out-dir ${coassembly_dir} ${CONFIG_FILEPATH} --until ${end_step} 2>&1 | tee ${coassembly_dir}/atlas_run_${log_code}.log

}


function find_atlas_binaries {
	# Description: finds location of binaries generated by ATLAS, to use for read mapping
	# GLOBAL Params: OUTPUT_DIR
	# Local params: none
	# Return: bbwrap_path (string); pileup_path (string); samtools_path (string)
	
	local coassembly_dir="${OUTPUT_DIR}/coassembly"
	
	bbwrap_path=$(find ${coassembly_dir}/.snakemake/conda -name "bbwrap.sh" | grep "/bin/bbwrap.sh")
	pileup_path=$(find ${coassembly_dir}/.snakemake/conda -name "pileup.sh" | grep "/bin/pileup.sh")
	samtools_path=$(find ${coassembly_dir}/.snakemake/conda -name "samtools" | grep "/bin/samtools")
	
	# TODO - add sanity check to make sure path was found and that only a single path was found. For now, just state what the paths were.
	echo "bbwrap_path: ${bbwrap_path}"
	echo "pileup_path: ${pileup_path}"
	echo "samtools_path: ${samtools_path}"
	
}


function read_map_to_coassemblies {
	# Description: iteratively maps read_mapping_samples to coassemblies like done within ATLAS
	# GLOBAL Params: OUTPUT_DIR; THREADS; coassembly_names (array); read_mapping_samples (array)
	# Local params: note that the following globals will be generated during running and will be used as inputs:  bbwrap_path (string); pileup_path (string); samtools_path (string)
	# Return: writes files/directories to disk

	echo "Read mapping iteratively to coassemblies..."

	local coassembly_dir="${OUTPUT_DIR}/coassembly"

	find_atlas_binaries # generates new variables with prefixes AND script name for binaries; see above

	# Manually add additional settings needed for scripts.
	# TODO - pull these settings (at least MEMORY) from the .yaml file!
	local MEMORY=65 #TODO
	local TMPDIR=${coassembly_dir}/tmp
	mkdir -p ${coassembly_dir}/tmp

	for i in $(seq 1 ${#coassembly_names[@]}); do
		# Set counter to be based on zero, not one
		local j=$((${i}-1))
		
		# Get name of coassembly
		local coassembly=${coassembly_names[${j}]}
		
		# Make relevant directories for storing output
		mkdir -p ${coassembly_dir}/${coassembly}/multi_mapping/unmapped_post_filter
		mkdir -p ${coassembly_dir}/${coassembly}/multi_mapping/logs
		mkdir -p ${coassembly_dir}/${coassembly}/multi_mapping/contig_stats
	
		# Temporarily change the internal fields separator (IFS) to parse comma separators. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
		local OFS="$IFS"
		IFS=,

		# Get names of individual samples provided for that coassembly name
		local mapping_sample_IDs=(${read_mapping_samples[${j}]})

		# Fix the IFS
		IFS="$OFS"
		
		echo "${coassembly_name}: ${mapping_sample_IDs[@]}"
		
		# Read map iteratively for each mapping ID
		for mapping in ${mapping_sample_IDs[@]}; do
			# TODO - pull more settings from .yaml file (these are FIXED right now)

			echo "rule align_reads_to_final_contigs (${mapping}):"
			local command=$(echo "${bbwrap_path} nodisk=t ref=${coassembly_dir}/${coassembly}/${coassembly}_contigs.fasta in1=${OUTPUT_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_R1.fastq.gz,${OUTPUT_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_se.fastq.gz in2=${OUTPUT_DIR}/${mapping}/sequence_quality_control/${mapping}_QC_R2.fastq.gz,null trimreaddescriptions=t outm=${coassembly_dir}/${coassembly}/multi_mapping/${mapping}.sam outu1=${coassembly_dir}/${coassembly}/multi_mapping/unmapped_post_filter/${mapping}_unmapped_R1.fastq.gz,${coassembly_dir}/${coassembly}/multi_mapping/unmapped_post_filter/${mapping}_unmapped_se.fastq.gz outu2=${coassembly_dir}/${coassembly}/multi_mapping/unmapped_post_filter/${mapping}_unmapped_R2.fastq.gz,null threads=${THREADS} pairlen=1000 pairedonly=t mdtag=t xstag=fs nmtag=t sam=1.3 local=t ambiguous=best secondary=t ssao=t maxsites=10 -Xmx${MEMORY}G 2> ${coassembly_dir}/${coassembly}/multi_mapping/logs/contig_coverage_stats_${mapping}.log")
			echo $command
			$command
			echo ""

			echo "rule pileup (${mapping}):"
			local command=$(echo "${pileup_path} ref=${coassembly_dir}/${coassembly}/${coassembly}_contigs.fasta in=${coassembly_dir}/${coassembly}/multi_mapping/${mapping}.sam threads=${THREADS} -Xmx${MEMORY}G covstats=${coassembly_dir}/${coassembly}/multi_mapping/contig_stats/postfilter_coverage_stats_${mapping}.txt hist=${coassembly_dir}/${coassembly}/multi_mapping/contig_stats/postfilter_coverage_histogram_${mapping}.txt basecov=${coassembly_dir}/${coassembly}/multi_mapping/contig_stats/postfilter_base_coverage_${mapping}.txt.gz concise=t physcov=t secondary=f bincov=${coassembly_dir}/${coassembly}/multi_mapping/contig_stats/postfilter_coverage_binned_${mapping}.txt 2>> ${coassembly_dir}/${coassembly}/multi_mapping/logs/contig_coverage_stats_${mapping}.log")
			echo $command
			$command
			echo ""

			echo "rule convert_sam_to_bam (${mapping}):"
			mkdir -p ${coassembly_dir}/tmp/${coassembly}/multi_mapping/alignment/${coassembly}_${mapping}_tmp # TODO - delete later?
			local command=$(echo "${samtools_path} view -@ ${THREADS} -bSh1 ${coassembly_dir}/${coassembly}/multi_mapping/${mapping}.sam | ${samtools_path} sort -m 1536M -@ ${THREADS} -T ${coassembly_dir}/tmp/${coassembly}/multi_mapping/alignment/${coassembly}_${mapping}_tmp -o ${coassembly_dir}/${coassembly}/multi_mapping/${mapping}.bam -O bam -")
			echo $command
			$command

			# TODO delete temp files? What is done in ATLAS?
			rm ${coassembly_dir}/${coassembly}/multi_mapping/${mapping}.sam
			echo ""

		done
	
	done

	echo ""

}


# TODO long term, find a better solution for this (e.g., conda install) rather than using docker, which requires an internet connection. Also use specific metabat version
 function look_for_metabat {
 	# Description: installs docker to be able to use metabat if not already installed
 	# TODO - consider installing metabat if not there; also consider generating metabat_path variable
 	# GLOBAL Params: OUTPUT_DIR; #TODO
 	# Local params: none
 	# Return: none

	# TODO CHECK if docker isn't installed!
	# TODO assumes root user!
if 

apt-get remove -y docker docker-engine docker.io
apt-get update
apt-get install -y apt-transport-https ca-certificates curl software-properties-common gnupg2
curl -fsSL https://download.docker.com/linux/$(. /etc/os-release; echo "$ID")/gpg | apt-key add -
add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/$(. /etc/os-release; echo "$ID") $(. /etc/os-release; echo "$VERSION" | cut -d "(" -f 2 | cut -d ")" -f 1) stable"
apt-get update
apt-get install -y docker-ce


 	
 }


function organize_new_bins {
	# Description: After binning with metabat, moves the generated bins into the coassembly ATLAS folder to appear as though they were generated by MaxBin, to work with the rest of the ATLAS pipeline.
	# GLOBAL Params: OUTPUT_DIR; coassembly_names (array); #TODO
	# Local params: none
	# Return: moves files on disk
	
}


function bin_coassemblies {
	# Description: runs metabat to bin coassemblies based on multiple read mapping
	# GLOBAL Params: OUTPUT_DIR; THREADS; coassembly_names (array); read_mapping_samples (array)
	# Local params: none
	# Return: writes files/directories to disk
	
	# # Check metabat is installed
	# look_for_metabat
	
	echo "Binning coassemblies..."
	local coassembly_dir="${OUTPUT_DIR}/coassembly"
	
	# Manually add additional settings needed for scripts.
	# TODO - pull these settings (at least MEMORY) from the .yaml file!
	local ${MINCONTIG}=1000
	local ${MIN_BIN_SIZE}=200000

	for i in $(seq 1 ${#coassembly_names[@]}); do
		# Set counter to be based on zero, not one
		local j=$((${i}-1))
		
		# Get name of coassembly
		local coassembly=${coassembly_names[${j}]}
		
		# Make relevant directories for storing output
		local bin_output_dir="${coassembly_dir}/${coassembly}/multi_mapping/genomic_bins"
		mkdir -p ${bin_output_dir}
	
		# Temporarily change the internal fields separator (IFS) to parse comma separators. See Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 407 and corresponding Github page README at https://github.com/vsbuffalo/bds-files/tree/master/chapter-12-pipelines (accessed Nov 19, 2017)
		local OFS="$IFS"
		IFS=,

		# Get names of individual samples provided for that coassembly name
		local mapping_sample_IDs=(${read_mapping_samples[${j}]})

		# Fix the IFS
		IFS="$OFS"
		
		echo "${coassembly_name}"
		
		# Build array of BAM file locations - iteratively add each filepath
		local bam_filepaths=($(echo ""))	
			
		for mapping in ${mapping_sample_IDs[@]}; do
			local bam_filepaths=($(echo "${bam_filepaths[@]} ${coassembly_dir}/${coassembly}/multi_mapping/${mapping}.bam"))
		done
		
		# Run metabat
		docker run metabat/metabat:latest runMetaBat.sh --minContig ${MINCONTIG} --numThreads ${THREADS} \
--maxP 95 --minS 60 --maxEdges 200 --minCV 1 --minCVSum 1 --minClsSize ${MIN_BIN_SIZE} --seed 0 --unbinned \
${coassembly_dir}/${coassembly}/${coassembly}_contigs.fasta ${bam_filepaths[@]} | tee ${coassembly_dir}/${coassembly}/multi_mapping/logs/genomic_binning.log
		
		# Reorganize output into ATLAS folder, then delete temporary folder
		organize_new_bins
	
	done

}


function make_atlas_dirs_round2 {
	# Description: makes (empty) files and directories with correct time stamps to prepare ATLAS to finish running
	# GLOBAL Params: OUTPUT_DIR; coassembly_names (array)
	# Local params: none
	# Return: writes (empty) files/directories to disk

}


function finish_atlas {
	# Description: runs standard ATLAS pipeline, starting from read mapping onward to the end
	# GLOBAL Params: OUTPUT_DIR; THREADS; CONFIG_FILEPATH; # TODO
	# Local params: none
	# Return: writes files/directories to disk

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
	run_atlas_assemble

	# Extension of ATLAS: mapping and binning
	read_map_to_coassemblies
	bin_coassemblies

	# Finish running ATLAS with the new bins
	make_atlas_dirs_round2
	finish_atlas

	end_time=$(date)	
	echo ""
	echo "Coassembly finished. Output saved in ${OUTPUT_DIR}/coassembly."
	# TODO - remove the comment below once this is addressed.
	echo "(Note that there will be many empty files/directories in place to force ATLAS to run normally. Don't be alarmed.)"
	echo "Started at ${start_time} and finished at ${end_time}."
	echo ""

	echo "$(basename $0): finished."
	echo ""

}

main

