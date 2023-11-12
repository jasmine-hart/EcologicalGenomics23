# EcoGen HW 2 

# Q2 - * How does the filtration of transcripts by read depth in DESeq2 affect gene expression results? To address this question, you would need to choose 2-3 ways to filter transcripts by read depth, one of which should be the one we ran in class (75% of samples have at least 15 reads), test for differential expression for each, then compare results, (how many transcripts included, what proportion differentially expressed, PCAs). You could focus on F0 samples for this.

setwd("/Users/jasminehart/Documents/GitHub/EcologicalGenomics23/Transcriptomics/results/")

# Load the package
library(WGCNA) # weighted gene coexpression
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

library(DESeq2)
library(ggplot2)

library(pheatmap)

library(tidyverse)

install.packages("remotes")
remotes::install_github("kevinblighe/CorLevelPlot")

library(CorLevelPlot) 
library(gridExtra)

library(Rmisc) 

library(ggplot2)

# 1. Import the counts matrix and metadata and filter using DESeq2

countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)
# 1 to say its agnostic and not tied to anything without time table

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds15 <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples
# unique short reads

dds1 <- dds[rowSums(counts(dds) >= 10) >= 28,]
nrow(dds1) 

# [1] 28196, that have at least 10 reads (counts) in 75% of the samples

dds2 <- dds[rowSums(counts(dds) >= 20) >= 28,]
nrow(dds2) 
#[1] 23295, that have at least 20 reads

dds25 <- dds[rowSums(counts(dds) >= 25) >= 28,]
nrow(dds25) 

# [1] 23295, that have at least 20 reads (counts) in 75% of the samples
# rarer allelles or LD regions or long regions

# Run the DESeq model to test for differential gene expression <10 and <20

dds1 <- DESeq(dds1)

dds15 <- DESeq(dds15)

dds2 <- DESeq(dds2)

dds25 <- DESeq(dds25)

# Check the quality of the data by sample clustering and visualization
# The goal of transformation "is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low."

library("pheatmap")
library("vsn")

# this gives log2(n + 1)
ntd <- normTransform(dds15)
meanSdPlot(assay(ntd))

# Variance stabilizing transformation
vsd <- vst(dds15, blind=FALSE)
meanSdPlot(assay(vsd))


sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

###############################################################
# PCA to visualize global gene expression patterns

# First normalize the data using variance stabilization
vsd <- vst(dds15, blind=FALSE)

data <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

###########  

dataF0 <- subset(data, generation == 'F0')

dataF0 <- ggplot(dataF0, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

dataF0


png("PCA_F0.png", res=300, height=5, width=5, units="in")

ggarrange(dataF0, nrow = 1, ncol=1)

dev.off()

#Check the quality of the data by sample clustering and visualization
# The goal of transformation "is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low."

library("pheatmap")
library("vsn")

# this gives log2(n + 1)
ntd2 <- normTransform(dds2)
meanSdPlot(assay(ntd2))

# Variance stabilizing transformation
vsd2 <- vst(dds2, blind=FALSE)
meanSdPlot(assay(vsd2))


sampleDists2 <- dist(t(assay(vsd2)))

# First normalize the data using variance stabilization
vsd2 <- vst(dds2, blind=FALSE)

data2 <- plotPCA(vsd2, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar2 <- round(100 * attr(data,"percentVar"))


plot <- ggplot(data2, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  coord_fixed()

###########  

dataF0_2 <- subset(data2, generation == 'F0')

F0_2 <- ggplot(dataF0_2, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) +
  ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

F0_2


png("PCA_F0_2.png", res=300, height=5, width=5, units="in")

ggarrange(F0_2, nrow = 1, ncol=1)

dev.off()

library("ggpubr")


# by generation and treatment

resAM_OWA <- results(dds2, name="treatment_OWA_vs_AM", alpha=0.05)

resAM_OWA <- resAM_OWA[order(resAM_OWA$padj),]
head(resAM_OWA)  

summary(resAM_OWA)


resAM_OW <- results(dds15, name="treatment_OW_vs_AM", alpha=0.05)

resAM_OW <- resAM_OW[order(resAM_OW$padj),]
head(resAM_OW)  

summary(resAM_OW)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK


table(gsg$goodGenes)
table(gsg$goodSamples)


# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 



# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x



pca.var <- pca$sdev^2

pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# 3. Normalization ----------------------------------------------------------------------

colData <- row.names(sample_metadata)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order



# perform variance stabilization
dds_norm1 <- vst(dds1)

dds_norm15 <-vst(dds15)

dds_norm2 <- vst(dds2)

dds_norm25 <- vst(dds25)

# dds_norm <- vst(normalized_counts)

# get normalized counts for 10 
norm.counts <- assay(dds_norm1) %>%
  t()

# get normalized counts for 15
norm.counts15 <- assay(dds_norm15) %>%
  t()

# get normalized counts for 20 <
norm.counts2 <- assay(dds_norm2) %>%
  t()

# get normalized counts for 25 <
norm.counts25 <- assay(dds_norm25) %>%
  t()


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

# Call the network topology analysis function; for 15<
sft15 <- pickSoftThreshold(norm.counts15,
                          powerVector = power,
                          networkType = "signed",
                          verbose = 5)


# Call the network topology analysis function; for 20<
sft2 <- pickSoftThreshold(norm.counts2,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

# Call the network topology analysis function; for 20<
sft25 <- pickSoftThreshold(norm.counts25,
                          powerVector = power,
                          networkType = "signed",
                          verbose = 5)

# signed means it cares about up or downregulated or correlated together, showing those interactions 

sft.data <- sft$fitIndices

sft.data15 <- sft15$fitIndices

sft.data2 <- sft2$fitIndices

sft.data25 <- sft25$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

library(ggplot2)
library(gridExtra)
grid.arrange(a1, a2, nrow = 2)

# for original filtering

a15 <- ggplot(sft.data15, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a25 <- ggplot(sft.data15, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a15, a25, nrow = 2)

# >20 read length filtering 

a12 <- ggplot(sft.data2, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a22 <- ggplot(sft.data2, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a12, a22, nrow = 2)

a125 <- ggplot(sft.data25, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a225 <- ggplot(sft.data25, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a125, a225, nrow = 2)
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

norm.counts15[] <- sapply(norm.counts15, as.numeric)

norm.counts2[] <- sapply(norm.counts2, as.numeric)

norm.counts25[] <- sapply(norm.counts25, as.numeric)

soft_power <- 6
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package


# this step also takes a few minutes; ideally your maxBlockSize is larger than your number of genes to run the memory-intensive network construction all at once.

# raised number of genes to 40,000 bc 26k and 30k was not enough
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 40000,
                          minModuleSize = 30, 
                          reassignThreshold=0,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

# this step also takes a few minutes; ideally your maxBlockSize is larger than your number of genes to run the memory-intensive network construction all at once.
bwnet15 <- blockwiseModules(norm.counts15,
                          maxBlockSize = 26000,
                          minModuleSize = 30, 
                          reassignThreshold=0,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

# this step also takes a few minutes; ideally your maxBlockSize is larger than your number of genes to run the memory-intensive network construction all at once.

bwnet2 <- blockwiseModules(norm.counts2,
                          maxBlockSize = 26000,
                          minModuleSize = 30, 
                          reassignThreshold=0,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

bwnet25 <- blockwiseModules(norm.counts25,
                           maxBlockSize = 26000,
                           minModuleSize = 40, 
                           reassignThreshold=0,
                           TOMType = "signed",
                           power = soft_power,
                           mergeCutHeight = 0.25,
                           numericLabels = F,
                           randomSeed = 1234,
                           verbose = 3)

# TOMtype (Topological Overlap Matrix type) parameter - unsigned - doesn't consider positive/negative co-expression
# signed - when you want to consider the direction of co-expression interaction, e.g., activating or inhibiting
# WGCNA often uses a dendrogram-based approach to identify modules. The choice of the 
# height cut in the dendrogram can determine the number of modules. Selecting a higher
# cut height results in fewer, larger modules, while a lower cut height leads to more, 
# smaller modules.

cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

head(module_eigengenes)

# og filtering
module_eigengenes15 <- bwnet15$MEs

head(module_eigengenes15)

# >20 filtering
module_eigengenes2 <- bwnet2$MEs

head(module_eigengenes2)

# >20 filtering
module_eigengenes25 <- bwnet25$MEs

head(module_eigengenes2)

# get number of genes for each module >10
table(bwnet$colors)

# get number of genes for each module >15
table(bwnet15$colors)

# get number of genes for each module >20
table(bwnet2$colors)

# Plot the dendrogram and the module colors before and after merging underneath >10
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# Plot the dendrogram and the module colors before and after merging underneath >15
plotDendroAndColors(bwnet15$dendrograms[[1]], cbind(bwnet15$unmergedColors, bwnet15$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = ">15 read Filtered Dendrogram")

# Plot the dendrogram and the module colors before and after merging underneath >15
plotDendroAndColors(bwnet2$dendrograms[[1]], cbind(bwnet2$unmergedColors, bwnet2$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = ">20 read Filtered Dendrogram")
# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module

# 6A. Relate modules to traits --------------------------------------------------
# module trait associations

traits <- sample_metadata[, c(5,8,11,14,17)]


# Define numbers of genes and samples >10
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Define numbers of genes and samples >15
nSamples15 <- nrow(norm.counts15)
nGenes15 <- ncol(norm.counts15)

# Define numbers of genes and samples >20
nSamples2 <- nrow(norm.counts2)
nGenes2 <- ncol(norm.counts2)

module.trait.corr15 <- cor(module_eigengenes15, traits, use = 'p')
module.trait.corr.pvals15 <- corPvalueStudent(module.trait.corr15, nSamples)

# p is for pearsons

# visualize module-trait association as a heatmap

heatmap.data15 <- merge(module_eigengenes15, traits, by = 'row.names')

head(heatmap.data15)

heatmap.data15 <- heatmap.data15 %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data15)

CorLevelPlot(heatmap.data15,
             x = names(heatmap.data15)[12:16],
             y = names(heatmap.data15)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))

CorLevelPlot(heatmap.data2,
             x = names(heatmap.data2)[12:16],
             y = names(heatmap.data2)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"))


module.gene.mapping15 <- as.data.frame(bwnet15$colors) # assigns module membership to each gene
module.gene.mapping15 %>% 
  filter(`bwnet15$colors` == 'yellow') %>% 
  rownames()

groups15 <- sample_metadata[,c(3,1)]
module_eigengene.metadata15 <- merge(groups15, heatmap.data15, by = 'row.names')

#Create a summary data frame of a particular module eigengene information
MEyellow_summary15 <- summarySE(module_eigengene.metadata15, measurevar="MEyellow", groupvars=c("Generation","treatment"))

#Plot a line interaction plot of a particular module eigengene
ggplot(MEyellow_summary15, aes(x=as.factor(Generation), y=MEyellow, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=MEyellow-se, ymax=MEyellow+se), width=.15) +
  geom_line(aes(color=treatment, group=treatment, linetype = treatment)) +
  scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  xlab("Generation") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6))

# 6B. Intramodular analysis: Identifying driver genes ---------------

# Get top hub genes (genes with highest connectivity in the network) in WGCNA 
hubs15  <-  chooseTopHubInEachModule(norm.counts15, bwnet15$colors, type = "signed", omitColors = "")
hubs15

#   black 
# "TRINITY_DN24998_c0_g1::TRINITY_DN24998_c0_g1_i3::g.49112::m.49112" 
# blue 
# "TRINITY_DN6125_c0_g1::TRINITY_DN6125_c0_g1_i1::g.23836::m.23836" 
# brown 
# "TRINITY_DN1784_c0_g1::TRINITY_DN1784_c0_g1_i3::g.9203::m.9203" 
# green 
# "TRINITY_DN239_c0_g1::TRINITY_DN239_c0_g1_i12::g.1858::m.1858" 
# grey 
# "TRINITY_DN13235_c0_g2::TRINITY_DN13235_c0_g2_i3::g.38482::m.38482" 
# magenta 
# "TRINITY_DN4230_c0_g2::TRINITY_DN4230_c0_g2_i1::g.18061::m.18061" 
# pink 
# "TRINITY_DN865_c0_g1::TRINITY_DN865_c0_g1_i27::g.5072::m.5072" 
# purple 
# "TRINITY_DN1336_c0_g1::TRINITY_DN1336_c0_g1_i4::g.7523::m.7523" 
# red 
# "TRINITY_DN3251_c0_g1::TRINITY_DN3251_c0_g1_i25::g.14716::m.14716" 
# turquoise 
# "TRINITY_DN2215_c0_g1::TRINITY_DN2215_c0_g1_i1::g.11114::m.11114" 
# yellow 
# "TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.3643

### Plot Individual genes  to check! ### 

d15 <-plotCounts(dds15, gene="TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434", intgroup = (c("treatment","Generation")), returnData=TRUE)
d_summary15 <- summarySE(d15, measurevar = "count", groupvars=c("Generation","treatment"))

ggplot(d_summary15, aes(x=as.factor(Generation), y=count, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=count-se, ymax=count+se), width=.15) +
  geom_line(aes(color=treatment, group=treatment, linetype = treatment)) +
  scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  xlab("Generation") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6))



# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure15 <- cor(module_eigengenes15, norm.counts15, use = 'p')
module.membership.measure.pvals15 <- corPvalueStudent(module.membership.measure15, nSamples)


module.membership.measure.pvals15[1:10,1:10]

# Make a heat map of gene expressions within modules.
# Use the norm.counts matrix, subset based on module membership
t_norm.counts15 <- norm.counts15 %>% t() %>% as.data.frame()

# Yellow module
purple_transcripts15 <- module.gene.mapping15 %>% 
  filter(`bwnet15$colors` == 'purple') %>% 
  rownames()

t_norm.counts_purple15 <- t_norm.counts15 %>% 
  filter(row.names(t_norm.counts15) %in% purple_transcripts15)
purple15 <- t_norm.counts_purple15 - rowMeans(t_norm.counts_purple15)
df15 <- as.data.frame(colData(dds15)[,c("Generation","treatment")])

library(ggplot2)

#blue to yellow color scheme
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "black", "yellow"))(paletteLength)
myBreaks <- c(seq(min(t_norm.counts_purple15), 0, length.out=ceiling(paletteLength/2) + 1),seq(max(t_norm.counts_purple15)/paletteLength, max(t_norm.counts_purple15), length.out=floor(paletteLength/2)))

# make it numeric
t_norm.counts_purple15 <- as.data.frame(lapply(t_norm.counts_purple15,
                                          function(x) as.numeric(as.character(x))))

pheatmap(t_norm.counts_purple15, color = myColor, breaks = myBreaks,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Purple")

# Make a heat map of gene expressions within modules.
# Use the norm.counts matrix, subset based on module membership
t_norm.counts2 <- norm.counts2 %>% t() %>% as.data.frame()


# visualize module-trait association as a heatmap

heatmap.data2 <- merge(module_eigengenes2, traits, by = 'row.names')

head(heatmap.data2)

heatmap.data2 <- heatmap.data2 %>% 
  column_to_rownames(var = 'Row.names')

names(heatmap.data15)

names(heatmap.data2)

CorLevelPlot(heatmap.data2,
             x = names(heatmap.data2)[12:16],
             y = names(heatmap.data2)[1:11],
             col = c("blue1", "skyblue", "white", "pink", "red"),
             main = ">20 Filtering Correlation Plot")


module.gene.mapping2 <- as.data.frame(bwnet2$colors) # assigns module membership to each gene
module.gene.mapping2 %>% 
  filter(`bwnet2$colors` == 'yellow') %>% 
  rownames()

groups2 <- sample_metadata[,c(3,1)]
module_eigengene.metadata2 <- merge(groups2, heatmap.data2, by = 'row.names')

#Create a summary data frame of a particular module eigengene information
MEyellow_summary2 <- summarySE(module_eigengene.metadata2, measurevar="MEyellow", groupvars=c("Generation","treatment"))

#Plot a line interaction plot of a particular module eigengene
ggplot(MEyellow_summary2, aes(x=as.factor(Generation), y=MEyellow, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=MEyellow-se, ymax=MEyellow+se), width=.15) +
  geom_line(aes(color=treatment, group=treatment, linetype = treatment)) +
  scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  xlab("Generation") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6))

# 6B. Intramodular analysis: Identifying driver genes ---------------

# Get top hub genes (genes with highest connectivity in the network) in WGCNA 
hubs2  <-  chooseTopHubInEachModule(norm.counts2, bwnet2$colors, type = "signed", omitColors = "")
hubs2

# black 
# "TRINITY_DN12327_c0_g1::TRINITY_DN12327_c0_g1_i2::g.37336::m.37336" 
# blue 
# "TRINITY_DN1076_c0_g1::TRINITY_DN1076_c0_g1_i1::g.6134::m.6134" 
# brown 
# "TRINITY_DN3654_c0_g1::TRINITY_DN3654_c0_g1_i11::g.16052::m.16052" 
# green 
# "TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434" 
# greenyellow 
# "TRINITY_DN24998_c0_g1::TRINITY_DN24998_c0_g1_i3::g.49112::m.49112" 
# grey 
# "TRINITY_DN637_c0_g1::TRINITY_DN637_c0_g1_i15::g.4143::m.4143" 
# magenta 
# "TRINITY_DN865_c0_g1::TRINITY_DN865_c0_g1_i27::g.5072::m.5072" 
# pink 
# "TRINITY_DN239_c0_g1::TRINITY_DN239_c0_g1_i12::g.1858::m.1858" 
# purple 
# "TRINITY_DN6821_c0_g1::TRINITY_DN6821_c0_g1_i1::g.25753::m.25753" 
# red 
# "TRINITY_DN3687_c0_g4::TRINITY_DN3687_c0_g4_i4::g.16337::m.16337" 
# turquoise 
# "TRINITY_DN6125_c0_g1::TRINITY_DN6125_c0_g1_i1::g.23836::m.23836" 
# yellow 
# "TRINITY_DN8661_c0_g1::TRINITY_DN8661_c0_g1_i14::g.30312::m.30312" 

### Plot Individual genes  to check! ### 

d2 <-plotCounts(dds2, gene="TRINITY_DN239_c0_g1::TRINITY_DN239_c0_g1_i12::g.1858::m.1858", intgroup = (c("treatment","Generation")), returnData=TRUE)
d_summary2 <- summarySE(d2, measurevar = "count", groupvars=c("Generation","treatment"))

plot <- ggplot(d_summary2, aes(x=as.factor(Generation), y=count, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=count-se, ymax=count+se), width=.15) +
  geom_line(aes(color=treatment, group=treatment, linetype = treatment)) +
  scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  xlab("Generation") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(panel.grid.minor.y = element_blank(), legend.position = "none", plot.margin = margin(0,6,0,6))



# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure2 <- cor(module_eigengenes2, norm.counts2, use = 'p')
module.membership.measure.pvals2 <- corPvalueStudent(module.membership.measure2, nSamples2)


module.membership.measure.pvals2[1:10,1:10]

# Make a heat map of gene expressions within modules.
# Use the norm.counts matrix, subset based on module membership
t_norm.counts2 <- norm.counts2 %>% t() %>% as.data.frame()

# Yellow module
purple_transcripts2 <- module.gene.mapping2 %>% 
  filter(`bwnet2$colors` == 'purple') %>% 
  rownames()

t_norm.counts_purple2 <- t_norm.counts2 %>% 
  filter(row.names(t_norm.counts2) %in% purple_transcripts2)
purple2 <- t_norm.counts_purple2 - rowMeans(t_norm.counts_purple2)
df2 <- as.data.frame(colData(dds2)[,c("Generation","treatment")])

library(ggplot2)

#blue to yellow color scheme
paletteLength <- 50
myColor2 <- colorRampPalette(c("dodgerblue", "black", "yellow"))(paletteLength)
myBreaks2 <- c(seq(min(t_norm.counts_purple2), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_purple2)/paletteLength, max(t_norm.counts_purple2), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_purple2, color = myColor, breaks = myBreaks,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = ">20 Filtering Heatmap (Purple)")

# make it numeric
t_norm.counts_purple2 <- as.data.frame(lapply(t_norm.counts_purple2,
                                               function(x) as.numeric(as.character(x))))

# Make a heat map of gene expressions within modules.
# Use the norm.counts matrix, subset based on module membership
t_norm.counts2 <- norm.counts2 %>% t() %>% as.data.frame()

# Yellow module
yellow_transcripts2 <- module.gene.mapping2 %>% 
  filter(`bwnet2$colors` == 'yellow') %>% 
  rownames()

t_norm.counts_yellow2 <- t_norm.counts2 %>% 
  filter(row.names(t_norm.counts2) %in% yellow_transcripts2)

t_norm.counts_yellow2 <- t_norm.counts_yellow2 - rowMeans(t_norm.counts_yellow2)
df2 <- as.data.frame(colData(dds2)[,c("Generation","Treatment")])

t_norm.counts_yellow22 <- as.numeric(t_norm.counts_yellow2)

#blue to yellow color scheme
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "black", "yellow"))(paletteLength)
myBreaks <- c(seq(min(t_norm.counts_yellow2), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_yellow2)/paletteLength, max(t_norm.counts_yellow2), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_yellow2, color = myColor, breaks = myBreaks,show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Yellow")


