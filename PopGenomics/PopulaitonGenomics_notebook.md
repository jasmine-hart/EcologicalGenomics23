# Population Genomics

## Author: Jasmine Hart

### Affiliation: University of Vermont, Plant and Soil Science

### E-mail contact: jhart12@uvm.edu

### Start Date: 09/11/2023

### End Date: TBD

### Project Descriptions: This notebook will document my workflow on the bioinformatics of the Pop genomics section of Ecologgical Genomics (fall 2023)

# Table of Contents:

-   [Entry 1: 2023-09-11](#id-section1)
-   [Entry 2: 2023-09-13](#id-section2)
-   [Entry 3: 2023-09-18](#id-section3)
-   [Entry 4: 2023-09-20](#id-section4)
-   [Entry 5: 2023-09-25](#id-section5)
-   [Entry 6: 2023-09-27](#id-section6)
-   [Entry 7: 2023-10-02](#id-section7)
-   [Entry 8: 2023-10-04](#id-section8)
-   [Entry 9: 2023-10-09](#id-section9)
-   [Entry 10: 2023-10-11](#id-section10)
-   [Entry 11: 2023-10-16](#id-section11)
-   [Entry 12: 2023-10-18](#id-section12)
-   [Entry 13: 2023-10-23](#id-section13)
-   [Entry 14: 2023-10-25](#id-section14)
-   [Entry 15: 2023-10-30](#id-section15)
-   [Entry 16: 2023-11-01](#id-section16)
-   [Entry 17: 2023-11-06](#id-section17)
-   [Entry 18: 2023-11-08](#id-section18)


------    
<div id='id-section1'/> 

### Entry 1: 2023-09-11.   

- We reviewed the red spruce study system and the exome capture data
- We discussed the structure of fastq files (dna sequence, plus the Q scores)
- Using the program fastqc, we analyzed the quality of the sequencing runs


------    
<div id='id-section2'/>   


### Entry 2: 2023-09-13.  

- After discussing ht FastQC results, we saw good quality sequence data for most of the read length
- The initial 5 bp were discarded bc they had a more variable base frequencies, and the very end had slightly lower Q-scored
- Based on this we set up an analysis to trim the reads using the `fastp` program
- We ran the bash script `fastp.sh` for this
- We looked at the html files pruned by `fastp` and compared pre and post trimming -- things looked good!
- We ended the day setting up our read mapping of the trimmed and cleaned reads using `bwa`


------    
<div id='id-section3'/>   


### Entry 3: 2023-09-18.

- After conducting the SAM (sequence alignment) files we processed our map files
- Using `sambamba` and `samtools`
- We converted the sam files to bam format
- Sorting the bam files by its read coordinates
- Mark and remove PCR duplicates
- The sorted, cleaned and index files were moved to a place for easy/quick look up


------    
<div id='id-section4'/>   


### Entry 4: 2023-09-20.

- Using `flagstat` we calculated how well the reads mapped to the reference
- `depth` (a samtools command) to estimate depth coverage (average reads/site)
- Using those two commands in a loop to sample our individual populations (I had group 2020, for example)
- `awk` was used to format the output
- We started working on the `ANGSD` code in vim
- Running the code to estimate genotype likelihoods for our populations.


------    
<div id='id-section5'/>   


### Entry 5: 2023-09-25.

- Our code finished running and they were input into our local directories ("saf" means site allele frequency)
- We created the `ANGSD_doTheta.sh`to estimate SFS (site frequency spectrum) and nucleotide diversity stats 


------    
<div id='id-section6'/>   


### Entry 6: 2023-09-27.

- Using ANGSD and nucleotide diversity stats to caluclate divergence between populations (Fst) (red spruce and the black spruce reference genomes and the populations thru the NE)
- PCA and Admixture were used to visualize the population structure
- `pcANGSD` was used and the results were moved with Filezilla to R
- Using R script, we visualized the PCA and admixture with four figures


------    
<div id='id-section7'/>   


### Entry 7: 2023-10-02.

- We used `pcANGSD` to scan Fst for outliers and to identify significant outliers
- Also identifying minor allele loci
- Exporting this info to R in order to visualize this in an easily digestible manner.
- Also explored plantgenie for signicant loci

------    
<div id='id-section8'/>   


### Entry 8: 2023-10-04.
- We ran GEA association with imported BioClim data to investigate
- We used the same `pcANGSD` results with the bioclim info
- Plotting it using R
- Then transferring it back to bash to identify significant loci.

------    
<div id='id-section9'/>   


### Entry 9: 2023-10-09.
- Orienting us to the data set.
- Background and important info explained, to set the stage for subsequent data analysis.


------    
<div id='id-section10'/>   


### Entry 10: 2023-10-11.
- Using fastp again to clean the raw reads.
- Transferring them to our computers.



------    
<div id='id-section11'/>   


### Entry 11: 2023-10-16.
- We used BUSCO to check the completeness of our assembly.
- Salmon was used to map our reads to our de novo transcriptome assembly.


------    
<div id='id-section12'/>   


### Entry 12: 2023-10-18.
- We transferred this data to our personal devices and switched to R.
- DeSEQ2 was used to begin analyzing gene expression
- We imported the data set and then ran some basic statistical tests already loaded into R.

------    
<div id='id-section13'/>   


### Entry 13: 2023-10-23.
- We started using th DESeq package to define experimental design and set read length parameters.
- Utilzing heat maps and SD to check the quality of the data clusters (`pheatmap` and `vsn`)
- Then assembling this data into a simple PCA determined by genreation and treatment.

------    
<div id='id-section14'/>   


### Entry 14: 2023-10-25.
- Filtered the metadata and reads matrix by importing into DESeq.
- Then we started analyzing these with `WGCNA`
- This is the stage where I played with filter settings for HW2 (higher and lower).
- We looked for outlier loci using to methods: hierarchical clustering and a PCA.
- Setting the soft power threshold to cluster genes in a biologically meaningful way.
- Then visualizing the pick power using scatter plots based on r^2 value.
- Running WGCNA to get a dendrogram that explores the network properties of genes.


------    
<div id='id-section15'/>   


### Entry 15: 2023-10-30.
- We used eigenegenes to determine the number of genes for each module or functional categorization.
- We then linked the modules and trait association.
- Using a heatmap to demonstrate correlations between traits and modules.
- Zooming in on gene expressions within modules of interest.


------    
<div id='id-section16'/>   


### Entry 16: 2023-11-01.
- Performed GO analysis on our DESeq results by generation and treatment.


------  


------    
<div id='id-section16'/>   


### Entry 17: 2023-11-06.
- I missed class this week so I will reporting the SV unit tasks in the next entry.




------    
<div id='id-section18'/>   


### Entry 18: 2023-11-08.
- We filtered the bcf files (sorted by chromosome) using `bcftools`
- Bcf files were converted to vcf files in the same step.
- `Local PCA` was then run to map the filtered chromosomes to the host genome.
- Regions (gene clusters/corners) of interest were further investigated using GO enrichment.
- Corners of the MDS corresponded to important genomic regions.


------  
