# atlas-extensions
Extensions to the ATLAS pipeline

Copyright Jackson M. Tsuji, 2018

**NOTE: These scripts are still in progress - early development only.

## Overview

[ATLAS](https://github.com/pnnl/atlas) is a promising, extendable pipeline for metagenomic analysis. This repo contains a collection of scripts to apply custom extensions on top of the existing ATLAS framework for personal analyses.

## Project 1: coassembly extension

This is a rough workaround to support coassembly and differential abundance binning using ATLAS release 1.0.22. The ATLAS team appears to be developing a full coassembly-like module that could replace this collection of scripts long-term.

The coassembly extension included in this repo only works for a subset of ATLAS settings. It requires that paired-end reads were used for metagenomic sequencing and hard-codes some of the settings for read mapping that would normally be set in the config (.yaml) file. Only use this extension with these caveats in mind.

It is recommended to use this extension within a Docker container (using the provided wrapper `enter-atlas-coassembly`) to avoid issues during install and usage.

**Dependencies**: [Docker](https://docs.docker.com/install/#supported-platforms)

**Requirements**: Requires that you have already run ATLAS release 1.0.22 normally on your metagenomes. See also the hard-coded settings above.

### Usage example
```
# Install enter-atlas-coassembly
cd [your_working_directory]
git clone https://github.com/jmtsuji/atlas-extensions.git
chmod 755 atlas-extensions/enter-atlas-coassembly
sudo mv atlas-extensions/enter-atlas-coassembly /usr/local/bin

# Enter the Docker container (first time entry could take several minutes)
[sudo] enter-atlas-coassembly [path_to_database_directory] [path_to_output_directory]
# See https://github.com/jmtsuji/atlas-wrapper for more detailed usage notes
# path_to_output_directory: must be the same as the output directory where the standard ATLAS run was performed

# Generate a config (.yaml) file for the run
co-assembly.sh /home/atlas/output /home/atlas/databases [CONFIG_FILEPATH.yaml] [COASSEMBLY_GUIDE_FILEPATH.tsv] [THREADS]
# CONFIG_FILEPATH.yaml: the path where you want the new coassembly .yaml file to be saved by co-assembly.sh
# COASSEMBLY_GUIDE_FILEPATH.tsv: tab-separated file indicating which samples from the standard ATLAS run should be coassembled and which samples should be read mapped for genome binning. See example file in this git repo under 'examples'.

# Modify the .yaml file to your liking, then perform the run
log_code=$(date '+%y%m%d_%H%M')
co-assembly.sh /home/atlas/output /home/atlas/databases [CONFIG_FILEPATH.yaml] [COASSEMBLY_GUIDE_FILEPATH.tsv] [THREADS] 2>&1 | tee /home/atlas/output/co-assembly_${log_code}.log
# CONFIG_FILEPATH.yaml: your modified version of the .yaml file output by the last step

# All done! Output will be saved to /home/atlas/output/coassembly

```

### More detailed usage notes
The concept of this extension is fairly simple. Script progression is as follows:
1. Combines QC-finalized reads from samples from the individual ATLAS run into input files for coassembly.
2. Runs ATLAS from the assembly step (megahit or spades, your choice) until the end of the pipeline. The snakefiles modules had to be slightly modified to omit QC summary rules because QC was not directly performed on the coasembly reads.
3. Iteratively read maps the read_mapping_samples onto the coassembled contigs using the same style of code as used in the standard ATLAS pipeline (hard-coded settings for now).
4. Performs differential abundance genome binning using [metabat2](https://bitbucket.org/berkeleylab/metabat). Replaces maxbin bins.
5. Re-runs the bin analysis steps of the ATLAS pipeline on the metabat2 bins. Integrates these into the final annotations.txt file. Also includes count stats of read mapping of indivdual samples onto the coassembly in the annotations.txt file. This finalized annotations file is saved as `[sample_ID}_annotations_multi_mapping.txt`.


