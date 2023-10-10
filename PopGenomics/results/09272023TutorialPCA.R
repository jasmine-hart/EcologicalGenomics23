setwd("/Users/jasminehart/Documents/GitHub/EcologicalGenomics23/PopGenomics/results/")

ggplot()
list.files()

## First, let's work on the genetic PCA:

COV <- as.matrix(read.table("allRS_poly.cov")) # read in the genetic covariance matrix

PCA <- eigen(COV) # extract the principal components from the COV matrix there are 95 PCs

## How much variance is explained by the first few PCs?

var <- round(PCA$values/sum(PCA$values),2)

var[1:2]

# A "screeplot" of the eigenvalues of the PCA:

barplot(var, 
        xlab="Eigenvalues of the PCA", 
        ylab="Proportion of variance explained")

## Bring in the bam.list file and extract the sample info:

names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

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

