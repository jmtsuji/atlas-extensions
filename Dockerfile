# Requires that pnnl/atlas:1.0.22 be built as: 'docker build -t pnnl/atlas:1.0.22 github.com/jmtsuji/atlas-wrapper.git'
FROM pnnl/atlas:1.0.22
LABEL maintainer="Jackson M. Tsuji <jackson.tsuji@uwaterloo.ca>"

# Install dependencies
RUN /bin/bash -c "conda install -y --name atlas_env -c r r-dplyr"
RUN git clone https://github.com/jmtsuji/atlas-extensions.git /home/atlas
RUN /bin/bash /home/atlas/atlas-extensions/setup/install_metabat2.sh /usr/local/bin
RUN mv /home/atlas/atlas-extensions/merge_atlas_multi_mapped_counts.R /usr/local/bin

# Add customized snakefile modules
RUN mv /home/atlas/atlas-extensions/setup/qc.snakefile /home/atlas/atlas-extensions/setup/assemble.snakefile /opt/conda/envs/atlas_env/lib/python3.6/site-packages/atlas/rules

ENTRYPOINT cd /home/atlas && \
	echo "Welcome to the ATLAS-coassembly docker container. Run 'atlas -h' to see run options. Type 'exit' to leave the container." && \
	/bin/bash
