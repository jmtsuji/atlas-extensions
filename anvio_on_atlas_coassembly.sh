# Running anvi'o from atlas-coassembly output
# Jackson M. Tsuji, Neufeld lab (Apr. 10, 2018)_
# All based off the Anvi'o metagenomics tutorial at http://merenlab.org/2016/06/22/anvio-tutorial-v2/ (etc.)
# Recommended to start from within a conda env

atlas_dir="/Hippodrome/jmtsuji/180123_ELA111314_atlas_r1.0.22_plus_full" # with "coassembly" dir inside
coassembly_sample_ID="CA-L227-2014"
output_dir=${atlas_dir}/post_analysis/04_anvio
threads=12

cogs_data_dir="/Hippodrome/anvio/databases"

mkdir -p ${output_dir}/01_import

# Test
# TODO - finish
if [ ! -d ${cogs_data_dir}/[database??] ]; then
	echo "ERROR: could not find COGs database at '${cogs_data_dir}/[database??]'. Exiting..."
	exit 1
fi

# 1. Prepare to import prokka info
cd ${output_dir}/01_import
wget https://raw.githubusercontent.com/karkman/gff_parser/master/gff_parser.py -O gff_parser.py
pip install gffutils

python gff_parser.py ${atlas_dir}/coassembly/${coassembly_sample_ID}/annotation/prokka/${coassembly_sample_ID}.gff \
				--gene-calls ${coassembly_sample_ID}_gene_calls.txt \
				--annotation ${coassembly_sample_ID}_gene_annot.txt

# 2. Generate contigs database
cd ${output_dir}
anvi-gen-contigs-database -f ${atlas_dir}/coassembly/${coassembly_sample_ID}_contigs.fasta \
				-o ${coassembly_sample_ID}_contigs.db -n ${coassembly_sample_ID}_contigs_db \
				--external-gene-calls ${output_dir}/01_import/${coassembly_sample_ID}_gene_calls.txt

## Example
# gene_callers_id	contig	start	stop	direction	partial	source	version
# 1	contig_01	1113	1677	f	0	program	v1.0
## "The statement above means that the index of the first nucleotide in any contig should be 0. In other words, we start counting from 0, not from 1."


# 3. Annotate with single-copy marker genes
cd ${output_dir}
anvi-run-hmms -c contigs.db --num-threads ${threads}
# anvi-run-hmms -c contigs.db --num-threads ${threads} --hmm-profile-dir [your_dir_for_custom_hmms]

anvi-display-contigs-stats contigs.db

# 4. Add functional info
# TODO - find a way to test whether or not setup is needed. Assumes already set up for now.
# anvi-setup-ncbi-cogs --cog-data-dir ${cogs_data_dir} --num-threads ${threads}

anvi-run-ncbi-cogs --cog-data-dir ${cogs_data_dir} --num-threads ${threads}



anvi-import-functions -c contigs.db -i gene_calls.txt

# Import functions - functional
gene_callers_id	source	accession	function	e_value
1	Pfam	PF01132	Elongation factor P (EF-P) OB domain	4e-23

As you can see,
Not every gene call has to be present in the matrix,
It is OK if there are multiple annotations from the same source for a given gene call,
It is OK if a give gene is annotated only by a single source.


anvi-import-taxonomy -c CONTIGS.db -i input_matrix.txt -p default_matrix

# Taxonomy
gene_callers_id	t_phylum	t_class	t_order	t_family	t_genus	t_species
1	 	 	 	 	 	Bacteroides fragilis
2	 	 	 	 	 	Bacteroides fragilis



anvi-profile -i SAMPLE-01.bam -c contigs.db --output-dir --sample-name

anvi-merge SAMPLE-01/PROFILE.db SAMPLE-02/PROFILE.db SAMPLE-03/PROFILE.db -o SAMPLES-MERGED -c contigs.db --skip-concoct-binning

# Import my bins
anvi-import-collection binning_results.txt -p SAMPLES-MERGED/PROFILE.db -c contigs.db --source "SOURCE_NAME"[metabat2] --contigs-mode --bins-info [example_bins_info_file.txt]

# external_binning_of_contigs.txt
204_10M_MERGED.PERFECT.gz.keep_contig_878	Bin_2
204_10M_MERGED.PERFECT.gz.keep_contig_6515	Bin_3
204_10M_MERGED.PERFECT.gz.keep_contig_1720	Bin_4
contig_1	Bin_4

# example_bins_info_file.txt
Bin_0	UNKNOWN_SOURCE	#ABCDEF
Bin_1	UNKNOWN_SOURCE	#FF2244
Bin_2	UNKNOWN_SOURCE	#2244FF
Bin_3	UNKNOWN_SOURCE	#22FF44
Bin_4	UNKNOWN_SOURCE	#44FFFF


anvi-interactive -p SAMPLES-MERGED/PROFILE.db -c contigs.db -C metabat2 --server-only -P 8080



