FROM continuumio/miniconda:4.4.10
LABEL maintainer="Jackson M. Tsuji <jackson.tsuji@uwaterloo.ca>"

# Update conda to prevent gcc error
RUN conda update -y conda python

# Install metabat2
RUN apt-get update
RUN git clone -b improved_docker https://github.com/jmtsuji/atlas-extensions.git /home/atlas/atlas-extensions
RUN /bin/bash /home/atlas/atlas-extensions/setup/install_metabat2.sh /usr/local /home/atlas/tmp

# Make conda env and install ATLAS and R dependencies
RUN conda create -y --name atlas_coassembly_env python=3.6
RUN conda install -y --name atlas_coassembly_env -c bioconda python=3.6 snakemake bbmap=37.17 click
RUN conda install -y --name atlas_env -c r r-plyr r-dplyr r-getopt
RUN /bin/bash -c "source activate atlas_env && pip install -U 'pnnl-atlas==1.0.22' && source deactivate"

# Install coassembly scripts
RUN mv /home/atlas/atlas-extensions/co-assembly.sh /home/atlas/atlas-extensions/merge_atlas_multi_mapped_counts.R /usr/local/bin

# Add customized snakefile modules
RUN mv /home/atlas/atlas-extensions/setup/qc.snakefile /home/atlas/atlas-extensions/setup/assemble.snakefile /opt/conda/envs/atlas_env/lib/python3.6/site-packages/atlas/rules

# Clean up
RUN rm -r /home/atlas/atlas-extensions

# Add lines to Snakefile to force bash shell usage if not already in ATLAS. Can remove in future version of ATLAS.
RUN if ! grep -q "shell\.executable" /opt/conda/envs/atlas_env/lib/python3.6/site-packages/atlas/Snakefile; then \
	sed -i '12i\shell.executable("/bin/bash")\nshell.prefix("set -o pipefail; ")' \
	/opt/conda/envs/atlas_env/lib/python3.6/site-packages/atlas/Snakefile; fi

# Add code to automatically start atlas_env environment when logging in
RUN echo "source activate atlas_coassembly_env" >> /root/.bashrc

RUN mkdir -p /home/atlas

ENTRYPOINT cd /home/atlas && \
	echo "Welcome to the ATLAS-coassembly docker container. Run 'co-assembly.sh' to see run options. Type 'exit' to leave the container." && \
	/bin/bash
