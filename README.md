# RNAseq_preprocessing
These are tutorials on a subset of tools available for processing raw RNAseq data. This is for two pipelines HISAT2_SAMtools_Stringtie_gffcompare_ballgown pipeline or HISAT2_SAMtools_Stringtie_PrepDEanalysis.py_DESeq2 pipeline. The first three steps are command line arguments and require installation of the following softwares.

`FASTQC`
`HISAT2`
`SAMTOOLS`
`Stringtie`
`gffcompare`


The description and tutorial of the tools was originally published here: https://www.nature.com/articles/nprot.2016.095
The sequenced samples used in tutorial was from original material PMID:24037378 or https://www.ncbi.nlm.nih.gov/pubmed/24037378 and is accessible.

The sample data used for the analysis in tutorial is a subset of samples used in the original study. And the subset of samples (12) are further filtered for reads that map onto human chromosome X, for ease of analysis and reducing computational demand. So the index files(ex. chrX_tran.1.ht2...), [which are used by HISAT2 for mapping] and reference annotated transcriptome (ex. chrX.gtf)[which is used by Stringtie for reference guided transcriptome assembly] are also of chromosome X. It is possible to expand the region covered, see below.

# Customizing the reference files required by aligners and assemblers
### It is possible to generate the index files yourself (computationally intensive!)
Although the reference genome is not used in this analysis it can be downloaded as follows
`curl -L "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" > hg38.fa.gz`

### It is possible to download the whole reference genes (refseq or refgene) from the UCSC genome browser

Start at: https://genome.ucsc.edu/cgi-bin/hgTables and then choose as follows
* clade: Mammal
* genome: Human
* assembly: Dec. 2013 (GRCh38/hg38)
* group: Genes and Gene Predictions
* track: NCBI RefSeq
* region: genome (radio button)
* table: UCSC RefSeq(refGene)
* output format: GTF-gene transfer format(limited)
* output file: GRCh38_hg38_ucsc_annotated.gtf
* Send output to: radio buttons unchecked.
* file type returned: plain text

** Finally press "get output" button.

https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=764817289_W7uU8Aj03hotPfShgOJGHr2RDr1M&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=refS$

### To get the complete list of samples with paired end reads
can ftp to the following website
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/

