#!/bin/bash
### Installing metabat2 in the miniconda container
# Jackson M. Tsuji (Neufeld lab PhD student)
# Rough script... version 1.0.0 (started Feb. 23rd, 2018)

# Basic script stuff (from Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 397):
set -e
set -u
set -o pipefail

script_version="1.0.1"

# If no input is provided, exit out and provide help
if [ $# == 0 ]
    then
    printf "$(basename $0): installs metabat on debian (jessie) to specified directory.\n"
    printf "Version: ${script_version}\n"
    printf "Contact Jackson M. Tsuji (jackson.tsuji@uwaterloo.ca) for bug reports or feature requests.\n\n"
    printf "Usage: $(basename $0) install_directory download_directory\n\n"
    printf "Usage details:\n"
    printf "1. install_directory: Path to the base directory where metabat binaries will be installed. The binaries will be saved in a subdirectory of the base directory called 'bin'. For example, you could choose '/usr/local'.\n"
    printf "2. download_directory: Path to the temporary directory where files involved in installation will be stored until the end of the script. Directory must not already exist and will be cleaned up at the end of the install.\n\n"
    printf "Install dependencies: apt; an internet connection.\n\n"
    printf "Final notes:\n"
    printf "* Installs dependency packages GLOBALLY via apt -- careful!"
    printf "* This script may also work for other linux distro's, but this has not been tested.\n\n"
    exit 1
fi
# Using printf: http://stackoverflow.com/a/8467449 (accessed Feb 21, 2017)
# Test for empty variable: Bioinformatics Data Skills Ch. 12 pg 403-404, and http://www.tldp.org/LDP/Bash-Beginners-Guide/html/sect_09_07.html and http://stackoverflow.com/a/2428006 (both accessed Feb 21, 2017)

# Set variables from user input:
INSTALL_DIR=$1
TEMP_DIR=$2

start_time=$(date)

echo "Running $(basename $0) on ${start_time}. Input: install_directory = '${INSTALL_DIR}'; download_directory = '${TEMP_DIR}'."

# Install dependencies
echo "Updating apt..."
apt-get update > /dev/null
echo "Installing dependencies..." # gcc binutils
apt-get install -y libboost-all-dev zlib1g-dev scons build-essential git curl libncurses5-dev
# conda update -y conda python
# conda install -y -c bioconda samtools # Seems like this isn't needed.

echo "Downloading metabat..."
metabat_dir=$(realpath ${TEMP_DIR})

if [ -d $metabat_dir ]; then
	echo "Specified temporary download folder '${metabat_dir}' already exists. Please delete before running this script. Job terminating."
	exit 1
fi

mkdir -p ${metabat_dir}
wget -P ${metabat_dir} https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz
cd ${metabat_dir}
tar -xzf master.tar.gz && rm master.tar.gz
cd berkeleylab-metabat-*

# Temporarily allow for unset variables in case BOOST_ROOT is null
echo "Installing and testing metabat2..."
set +u
scons install PREFIX=${INSTALL_DIR} install
set -u

cd $INSTALL_DIR
rm -r ${metabat_dir}

end_time=$(date)

echo ""
echo ""
echo "Done installing metabat2 to '${INSTALL_DIR}'"
echo "Started at ${start_time} and finished at ${end_time}."
echo ""

