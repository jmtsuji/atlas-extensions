#!/bin/bash
### Installing metabat2 in the miniconda container
# Jackson M. Tsuji (Neufeld lab PhD student)
# Rough script... version 1.0.0 (started Feb. 23rd, 2018)

# Basic script stuff (from Vince Buffalo's "Bioinformatics Data Skills" (1st Ed.) chapter 12, pg 397):
set -e
set -u
set -o pipefail

echo "Installing metabat2..."

starting_dir=$PWD

# Install dependencies
apt-get update && apt-get install -y libboost-all-dev zlib1g-dev scons build-essential git curl libncurses5-dev # gcc binutils
# Assuming python already there...
conda install -y -c bioconda samtools # maybe this isn't needed?

metabat_dir="/home/atlas/output/coassembly/.metabat"
mkdir -p ${metabat_dir}
wget -P ${metabat_dir} https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz
cd ${metabat_dir}
tar -xvzf ${metabat_dir}/master.tar.gz && rm ${metabat_dir}/master.tar.gz
cd ${metabat_dir}/berkeleylab-metabat-*
scons install PREFIX=$HOME [BOOST_ROOT=$BOOST_ROOT]

cd $starting_dir

echo "Done installing to ${metabat_dir}/berkeleylab-metabat-*"

