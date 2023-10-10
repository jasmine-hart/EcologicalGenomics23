library(RcppCNPy)# for reading python numpy (.npy) files

setwd("~/Documents/GitHub/EcologicalGenomics23/PopGenomics/results/")

list.files()

### read in selection statistics (these are chi^2 distributed)

s<-npyLoad("allRS_poly.selection.npy")

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,1))
names(pval) = c("p_PC1","p_PC2")

## read positions
p <- read.table("allRS_poly_mafs.sites",sep="\t",header=T, stringsAsFactors=T)
dim(p)

p_filtered = p[which(p$kept_sites==1),]

dim(p_filtered)

# get all the outliers with p-values below some cutoff
cutoff=1e-3   

outliers_PC1 <- p_filtered[which(pval$p_PC1<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff?
dim(outliers_PC1)[1]


# write them out to a file
write.table(outliers_PC1,
            "allRS_poly_outliers_PC1.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F)

COV <- as.matrix(read.table("allRS_poly.cov"))

PCA <- eigen(COV)

data=as.data.frame(PCA$vectors)
data=data[,c(1:2)] # the second number here is the number of PC axes you want to keep

write.table(data,
            "allRS_poly_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)

setwd("~/Documents/GitHub/EcologicalGenomics23/PopGenomics/results/")

bio <- getData("worldclim",var="bio",res=10)

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)

#The chunk below refers to your bamlist file that you transferred during last week's PCA/admixture analysis.  It should be the same one you want to use here -- if your sample list for analysis changes in the future, you'll need a different bamlist!

names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

angsd_coords <- merge(pops, coords, by.x="Ind", by.y="Tree")

points <- SpatialPoints(angsd_coords[c("Longitude","Latitude")])

clim <- extract(bio,points)

angsd_coords_clim <- cbind.data.frame(angsd_coords,clim)
str(angsd_coords_clim)

# Make the climate PCA:

clim_PCA = PCA(angsd_coords_clim[,15:33], graph=T)

# Get a screeplot of cliamte PCA eigenvalues

fviz_eig(clim_PCA)

# What is the climate PCA space our red spruce pops occupy?

fviz_pca_biplot(clim_PCA, 
                geom.ind="point",
                col.ind = angsd_coords_clim$Latitude, 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                title="Climate PCA (Bioclim)",
                legend.title="Latitude")

# Which variables show the strongest correlation on the first 2 climate PC axes?

dimdesc(clim_PCA)[1:2]

# Replace "XX" with your bio variable most significant on climate PC1:

write.table(scale(angsd_coords_clim["bio12"]),
            "allRS_bio12.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)


# Replace "YY" with your bio variable most significant on climate PC2:  

write.table(scale(angsd_coords_clim["bio10"]),
            "allRS_bio10.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)


COV3 <- as.matrix(read.table("allRS_poly3.cov")) # read in the genetic covariance matrix

PCA3 <- eigen(COV3) # extract the principal components from the COV matrix there are 95 PCs

## How much variance is explained by the first few PCs?

var3 <- round(PCA$values/sum(PCA$values),2)

var3[1:3]

# A "screeplot" of the eigenvalues of the PCA:

barplot(var3, 
        xlab="Eigenvalues of the PCA", 
        ylab="Proportion of variance explained")

## Bring in the bam.list file and extract the sample info:

names3 <- read.table("allRS_bam3.list")
names3 <- unlist(strsplit(basename(as.character(names[,0])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops3 <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names3(pops) = c("Ind", "Pop", "Row", "Col")

## A quick and humble PCA plot:

plot(PCA$vectors[,1:2],
     col=as.factor(pops[,2]),
     xlab="PC1",ylab="PC2", 
     main="Genetic PCA")

## A more beautiful PCA plot using ggplot :)

data=as.data.frame(PCA$vectors)
data=data[,c(1:3)]
data= cbind(data, pops)

cols=c("#377eB8","#EE9B00","#0A9396","#94D2BD","#FFCB69","#005f73","#E26D5C","#AE2012", "#6d597a", "#7EA16B","#d4e09b", "gray70")

ggscatter(data, x = "V1", y = "V2",
          color = "Pop",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))




## Next, we can look at the admixture clustering:

# import the ancestry scores (these are the .Q files)

q <- read.table("allRS_poly.admix.2.Q", sep=" ", header=F)

K=dim(q)[3] #Find the level of K modeled

## order according to population code
ord<-order(pops[,1])

# make the plot:
barplot(t(q)[,ord],
        col=cols[1:K],
        space=0,border=NA,
        xlab="Populations",ylab="Admixture proportions",
        main=paste0("Red spruce K=",K))
text(tapply(1:nrow(pops),pops[ord,2],mean),-0.05,unique(pops[ord,2]),xpd=T)
abline(v=cumsum(sapply(unique(pops[ord,3]),function(x){sum(pops[ord,3]==x)})),col=1:2,lwd=1.2)


# k=3 attempt
cols=c("#377eB8","#EE9B00","#0A9396","#94D2BD","#FFCB69","#005f73","#E26D5C","#AE2012", "#6d597a", "#7EA16B","#d4e09b", "gray70")
q3 <- read.table("allRS_poly.admix.3.Q", sep=" ",header=F)

K3=dim(q3)[2] #find the level of k modeled

## order according to pop code
ord<-order(pops[,2])

#make plot:
barplot(t(q3)[,ord],
        col=cols[1:K3],
        space=0,border=NA,
        xlab="Populations",ylab="Admixture proportions",
        main=paste0("Red Spruce K=",K3))
text(tapply(1:nrow(pops),pops[ord,2],mean),-0.05,unique(pops[ord,2]),xpd=T)
abline(v=cumsum(sapply(unique(pops[ord,3]),function(x){sum(pops[ord,3]==x)})),col=1:2,lwd=1.2)

# K3 Outlier loci

library(RcppCNPy)# for reading python numpy (.npy) files

setwd("~/Documents/GitHub/EcologicalGenomics23/PopGenomics/results/")

list.files()

### read in selection statistics (these are chi^2 distributed)

s3<-npyLoad("allRS_poly.selection3.npy")

# convert test statistic to p-value for k3
pval3 <- as.data.frame(1-pchisq(s3,1))
names(pval3) = c("p_PC1","p_PC2")

## read positions
p3 <- read.table("allRS_poly_mafs3.sites",sep="\t",header=T, stringsAsFactors=T, fill=T)
dim(p3)

p3_filtered = p3[which(p3$kept_sites==1),]
dim(p3_filtered)

# get all the outliers with p-values below some cutoff
cutoff=1e-3   

outliers_PC13 <- p3_filtered[which(pval3$p_PC1<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff?
dim(outliers_PC13)[1]

# How many sites got filtered out when testing for selection? Why?

## make manhattan plot
plot(-log10(pval3$p_PC1),
     col=p3_filtered$chromo,
     xlab="Position",
     ylab="-log10(p-value)",
     main="Selection outliers: pcANGSD e=3 (K3)")

# We can zoom in if there's something interesting near a position...

plot(-log10(pval3$p_PC1[2e05:2.01e05]),
     col=p3_filtered$chromo, 
     xlab="Position", 
     ylab="-log10(p-value)", 
     main="Selection outliers: pcANGSD e=3 (K3)")

# get the contig with the lowest p-value for selection
sel_contig3 <- p3_filtered[which(pval3==min(pval3$p_PC1)),c("chromo","position")]
sel_contig3

# get all the outliers with p-values below some cutoff
cutoff=1e-3   # equals a 1 in 5,000 probability
outlier_contigs3 <- p3_filtered[which(pval3<cutoff),c("chromo","position")]
outlier_contigs3

# how many outlier loci < the cutoff?
dim(outlier_contigs3)[1]

# how many unique contigs harbor outlier loci?
length(unique(outlier_contigs3$chromo))

#output to bash
write.table(unique(outlier_contigs3$chromo),
            "allRS_poly3_PC1_outlier_contigs.txt", 
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

# write them out to a file
write.table(outliers_PC13,
            "allRS_poly3_outliers_PC1.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F)

COV3 <- as.matrix(read.table("allRS_poly3.cov"))

PCA3 <- eigen(COV3)

data3=as.data.frame(PCA3$vectors)
data3=data3[,c(1:2)] # the second number here is the number of PC axes you want to keep

write.table(data3,
            "allRS_poly4_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)

library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(geodata)

setwd("~/Documents/GitHub/EcologicalGenomics23/PopGenomics/results/")

bio <- getData("worldclim",var="bio",res=10)

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)

#The chunk below refers to your bamlist file that you transferred during last week's PCA/admixture analysis.  It should be the same one you want to use here -- if your sample list for analysis changes in the future, you'll need a different bamlist!

names3 <- read.table("allRS_bam.list3")
names3 <- unlist(strsplit(basename(as.character(names3[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")

pops3 <- data.frame(names3[1:95], do.call(rbind, split[1:95]))
names(pops3) = c("Ind", "Pop", "Row", "Col")

angsd_coords <- merge(pops3, coords, by.x="Ind", by.y="Tree")

points3 <- SpatialPoints(angsd_coords[c("Longitude","Latitude")])

clim <- extract(bio,points3)

angsd_coords_clim <- cbind.data.frame(angsd_coords,clim)
str(angsd_coords_clim)

# Make the climate PCA:

clim_PCA = PCA(angsd_coords_clim[,15:33], graph=T)

# Get a screeplot of cliamte PCA eigenvalues

fviz_eig(clim_PCA)

# What is the climate PCA space our red spruce pops occupy?

fviz_pca_biplot(clim_PCA, 
                geom.ind="point",
                col.ind = angsd_coords_clim$Latitude, 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                title="Climate PCA (Bioclim)",
                legend.title="Latitude")

# Which variables show the strongest correlation on the first 2 climate PC axes?

dimdesc(clim_PCA)[1:2]

# Replace "XX" with your bio variable most significant on climate PC1:

write.table(scale(angsd_coords_clim["bio12"]),
            "allRS3_bio12.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)


# Replace "YY" with your bio variable most significant on climate PC2:  

write.table(scale(angsd_coords_clim["bio10"]),
            "allRS3_bio10.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)
