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
