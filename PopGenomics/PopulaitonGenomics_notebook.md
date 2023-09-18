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

------    
<div id='id-section1'/> 

### Entry 1: 2023-09-11.   

- We reviewed the red spruce study systen and the exome capture data
- We discussed the structure of fastq files (dna sequene, pkus the Q scores)
- Using the program fastqc, we analyzed the quality of the sequencing runs


------    
<div id='id-section2'/>   


### Entry 2: 2023-09-13.  

- After discussing ht FastQC results, we saw good quality sequence data for most of the read length
- The initial 5 bp were discarded bc they had a more variable base frequencies, and the very end had slightly lower Q-scored
- Based on this we set up an analysis to trim th ereads using the `fastp` program
- We ran the bash script `fastp.sh` for this
- We looked at the html files prud ed by `fastp` and compared pre and post trimming -- things looked good!
- We ended the day setting up our read mapping of the trimmed and cleaned reads using `bwa`


------    
<div id='id-section3'/>   


### Entry 3: 2023-XX-XX.