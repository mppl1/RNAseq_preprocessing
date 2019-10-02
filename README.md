# RNAseq_preprocessing
These are tutorials on a subset of tools available for processing raw RNAseq data. This is for two pipelines HISAT2_SAMtools_Stringtie_gffcompare_ballgown pipeline or HISAT2_SAMtools_Stringtie_PrepDEanalysis.py_DESeq2 pipeline. The first few steps are command line arguments and require installation of the following softwares. See Reference section for details on each software.

`FASTQC`
`HISAT2`
`SAMTOOLS`
`Stringtie`
`gffcompare`


The description and tutorial of the tools was originally published here: https://www.nature.com/articles/nprot.2016.095
The sequenced samples used in tutorial was from original material [Lappalainen et al. 2013 Nature](http://dx.doi.org/10.1038/nature12531) and is accessible[PMID:24037378](https://www.ncbi.nlm.nih.gov/pubmed/24037378). A useful portal to explore the dataset from original study is provided at https://www.ebi.ac.uk/Tools/geuvadis-das/. The dataset consists of lymphoblastic cell lines derived from 465 individuals. The individuals are either male or female and from one of five different populations: CEPH (CEU--Utah residents of Western and Northern European ancestry), Finns (FIN--Finland), British (GBR), Toscani (TSI--Italy) and Yoruba (YRI--Nigeria).

The sample data provided by JHU tools (HISAT2 and Stringtie) used for the analysis in tutorial(https://www.nature.com/articles/nprot.2016.095) is a subset of samples used in the original study. And the subset of samples (12) are further filtered for reads that map onto human chromosome X, for ease of analysis and reducing computational demand. So the index files provided by HISAT2 or Stringtie during installation(ex. chrX_tran.1.ht2...,which are used by HISAT2 for mapping) and reference annotated transcriptome provided by Stringtie during installation (ex. chrX.gtf, which is used by Stringtie for reference guided transcriptome assembly) are also of chromosome X(note: typically, males have 1 X chromosome). For details of installing the bioinformatic tools and running the example dataset, please refer to Dave Tang's excellent [tutorial](https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/). In this github project are presented options to expand the region covered, see below. And also a way to automate the pipeline (from aligning the reads, to sorting and conversion to bam file type, transcriptome assembly), a bash script is provided that can be modified.

# Customizing the reference files required by aligners and assemblers and expanded sample data
### To download the expanded dataset, several options are presented
**To get the complete list of samples with paired end RNA reads

* can ftp to the following website and perform a recursive download(limited by network capacity both locally and host's)
 `wget -cr ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/`
 
  **or individual samples for example: 
  `wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188021/ERR188021_1.fastq.gz`	

* At European Nucleotide Archive (ENA) search for ERP001942
https://www.ebi.ac.uk/ena/data/view/PRJEB3366
http://www.ebi.ac.uk/ena/data/view/ERR188021-ERR188482

* Or accessed at ArrayExpress with accession number E-GEUV-1 (There is also an miRNA and miRNA&RNA bundle)
https://www.ebi.ac.uk/arrayexpress/

** To download the phenodata or metadata on samples visit: 
https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/samples/?full=true
* or at command line
`curl -L "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt">geuvadis_complete_phenodata.txt`
**To process the downloaded datatable into a subset with just 3 columns of interest
 * `wc -l geuvadis_complete_phenodata.txt` 
 * `awk -F "\t" '{print $NF}' geuvadis_complete_phenodata.txt | wc -l`
 * `awk -F "\t" '{print $34 " " $7 " " $36 "  "}' geuvadis_complete_phenodata.txt >geuvadis_3columns.txt`
 * `awk -F "\t" '{print $1 " " $2 " " $3}' geuvadis_3columns.txt`
 * `sed -i 's/^Com.*$/ids " " sex " " population/' geuvadis_3columns.txt` 
 * or `sed -i 's/^Com.*$/ids \t sex \t population/' geuvadis_3columns.txt`

To be consistent with Dave Tang's or JHU tutorial the relevant fields in the table for our purposes are "Characteristics[sex]"==gender of individual, "Characteristics[ancestry category] or Factor Value[ancestry category]"==geographic name of population, "Comment[ENA_RUN]"==sample names, "Comment[FASTQ_URI]"==sample file name and ftp path for download
And thus you can construct a new data table with just the relevant information as described above using bash or Python and R below

** To generate the phenodata after downloading using Python or R
  * python approach
  
    * import pandas as pd
    * geuvadis=pd.read_csv("geuvadis_complete_phenodata.txt", sep="\t", header=0)
    * geuvadis.shape
    * **The above should output: (924, 37)
    * geuvadis.iloc[0,:]
    * geuvadis.columns
    * geuvadis.iloc[0,33] sample name
    * geuvadis.iloc[0,35] The geo
    * geuvadis.iloc[0,6] gender
  
    **subset the dataframe/datatable to contain only relevant columns/variables
    * geuvadis_phenodata=geuvadis.iloc[:, [33,6,35]]
  
    **rename the column names to be consistent
    * geuvadis_phenodata.columns=['ids', 'sex', 'population']
  
  **one can also further subset the phenodata to just two populations and abbreviate the geographic population label
  
    * to write out the data frame as a csv file
    * export_csv = df.to_csv (r'geuvadis_phenodata.csv', sep="," , index = None, header=True) 

  
  * R approach
  
    * geuvadisdf<-read_csv("geuvadis_complete_phenodata.txt", sep="\t", header=T)
    * dim(geuvadisdf)
    * geuvadisdf<-geuvadisdf[,c(34,7,36]  
    * names(geuvadisdf)<-c('ids', 'sex', 'population')
    * write.table(geuvadisdf, "geuvadisdf_phenodata.csv", sep=",")

### It is possible to generate the index files yourself (computationally intensive! 160GB of RAM for whole genome) using HISAT2
 ** Another option is to use the aligner: STAR

HISAT2 will use the reference genome for generating index files
Although the reference genome is not used in this tutorial analysis it can be downloaded as follows
`curl -L "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" > hg38.fa.gz`

** However, JHU has these indexes available ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz

### Once index files are generated or available, then the samples can be aligned or mapped using the index files
  **Here is an example alignment command with HISAT2 using the JHU sample dataset(descibed above)
  
`hisat2 -p 4 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188428_chrX_1.fastq.gz -2 chrX_data/samples/ERR188428_chrX_2.fastq.gz -S ERR188428_chrX.sam`

** [options]*
* -p <num_cores>: # of processor core threads
* --dta: downstream transcriptome assembly. Ex. for Stringtie etc.
** Main arguments
* -x:basename of the index for the reference genome, which excludes the ".1.ht2" etc.
* -1: comma-separated list of files containing mate 1s
* -2: comma-separated list of files containing mate 2s
* -S: file to write SAM alignments to. If not specified it is sent to stdout or console


#### It is possible to download the whole reference genes (refseq or refgene) from the UCSC genome browser

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
** Note that this reference genome was publicly available after the lymphoblastic cell line study was published

https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=764817289_W7uU8Aj03hotPfShgOJGHr2RDr1M&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=refS$

### To sort and convert the sam files into bam files usig samtools
`samtools sort -@ 4 -o ERR188428_chrX.bam ERR188428_chrX.sam`
**explanations of options arguments
* -@ 4: is the key value pair of CPU threads to use
* -o: is to specify the output filename
* last is the name of input file

### To perform transcriptome assembly with Stringtie
** Once samples are mapped or aligned onto a reference genome, then Stringtie can be used to assemble the transcriptome.
An example assembly command with Stringtie using the JHU sample dataset(descibed above)

`stringtie map/ERR188428_chrX.bam -l ERR188428 -p 4 -e -G chrX_data/genes/chrX.gtf -o assembly/ERR188428_chrX.gtf`

**explanations of options arguments**
* path to input file
* -l:prefix for output transcripts 
* -p: number of CPU threads to use 
* -G: Reference annotation to guide assemble in either GTF/GFF3 files
* -o: output file path and name
* -e: this option is needed for DESEQ2. read bundles with no match to reference will be skipped http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#deseq

**and subsequently to merge the transcripts (requires a mergelist text file)
`stringtie --merge -p 4 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt`

### To summarize the assembled transcript in comparison to reference genes
`gffcompare -r chrX_data/genes/chrX.gtf -G -o merged stringtie_merged.gtf`

# Installation of softwares in Linux Ubuntu 18.04LTS
** For convenience in using the binaries it is desirable, after installation, to link them to /usr/local/bin(this varies with OS you are using) so that they can be invoked without specifing their full path.
`ln -s /<full>/<path>/<to>/<file> /usr/local/bin`

### FASTQC (precompiled binary). It is useful for quality control of RNAseq of samples
https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
Downloaded FastQC v0.11.8 (Win/Linux zip file)
`unzip fastqc_v0.11.8.zip`
**make binary executable
`chmod 755 fastq`
** To launch GUI
`./fastqc`
** To link it to /usr/local/bin
`ln -s ~/Desktop/FastQC/fastqc /usr/local/bin `

### HISAT2 (precompiled binary)
`wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip`
`unzip hisat2-2.1.0-Linux_x86_64.zip`
**To link it to /usr/local/bin
`ln -s ~/Desktop/HISAT2/hisat2-2.1.0/* /usr/local/bin`

### Samtools (from source file, requires build. Alternatively, there is debian package but an earlier version)
`wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2`
`tar -tvf samtools-1.9.tar.bz2`
`cd samtools-1.9.tar.bz2`
`./configure` **another option is `./configure --prefix=/where/to/install`
`make`
** it took 40seconds
`sudo make install`
** To uninstall if needed
`sudo make uninstall`
** to check the install steps
`make -n install`
**Linking to /usr/local/bin is unnecessary as the installation procedure moves the necessary binaries to appropriate folder 

### Stringtie (precompiled binary)
`wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.0.1.Linux_x86_64.tar.gz`
`tar -xzvf stringtie-2.0.1.Linux_x86_64.tar.gz`
`cd stringtie-2.0.1.Linux_x86_64`
`ln -s ~/Desktop/stringtie-2.0.1.Linux_x86_64/stringtie /usr/local/bin`

### GFFCOMPARE
`git clone https://github.com/gpertea/gclib`
`git clone https://github.com/gpertea/gffcompare`
`cd gffcompare`
`make release`
`ln -s ~/Desktop/gffcompare/* /usr/local/bin`
** may want to remove the unnecessary links later or initially select specific files to link instead of whole folder

**Alternatively, a precompiled version can be downloaded
`wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.11.4.Linux_x86_64.tar.gz`
`tar xzvf gffcompare-0.11.4.Linux_x86_64.tar.gz`
`cd gffcompare-0.11.4.Linux_x86_64`



### PrepDE.py (JHU script to prepare a count matrix of Stringtie output for DESeq2 or EdgeR)
wget http://ccb.jhu.edu/software/stringtie/dl/prepDE.py
chmod 755 prepDE.py
**usage
`python2 prepDE.py -i sample_lst.txt`


# References
* Lappalainen T et al "Transcriptome and genome sequencing uncovers functional variation in humans". Nature(501):506–511 (2013.(http://dx.doi.org/10.1038/nature12531)
* Pertea M et al "Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown". Nature Protocols (11):1650–1667 (2016). https://www.nature.com/articles/nprot.2016.095
* https://www.ebi.ac.uk/Tools/geuvadis-das/
* https://www.internationalgenome.org/
* https://en.wikipedia.org/wiki/1000_Genomes_Project
* CEPH https://www.coriell.org/1/NIGMS/Collections/CEPH-Resources. https://en.wikipedia.org/wiki/Fondation_Jean_Dausset-CEPH
* Dave Tang's tutorial blog https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/
* FASTQC (version: 0.11.8) https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
* HISAT2(version: 2.1.0) https://ccb.jhu.edu/software/hisat2/index.shtml
* Stringtie (version: 2.0.1) http://ccb.jhu.edu/software/stringtie/index.shtml 
* Gffcompare (version: ) http://ccb.jhu.edu/software/stringtie/gff.shtml#gffcompare or https://github.com/gpertea/gffcompare

