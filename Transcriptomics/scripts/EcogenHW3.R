# Ecological genomics HW 3


title: "Local PCA results"
date: "`r date()`"
---
  
  ```{r setup, include=FALSE}
library(lostruct)
library(colorspace)
library(jsonlite)
library(RColorBrewer)
fig.dim <- 4
knitr::opts_chunk$set(fig.width=2*fig.dim,fig.height=fig.dim,fig.align='center')
# set do.pdfs to TRUE to output PDFs of figures as well
if (!exists("do.pdfs")) { do.pdfs <- TRUE }
```
```{r plot_setup, include=FALSE}
layout_heights <- function (k,dl=0,ncol=1) {
  # to set up layout without 'dl' lines between plots
  # use like layout(1:5,heights=layout_heights(5))
  if (k==1) return(1)
  layout(matrix(seq_len(k*ncol),ncol=ncol))  # this changes par("csi")
  ds <- dl*par("lheight")*par("csi")
  eps=par("mai")[c(1,3)]
  dh=(par("din")[2]-sum(eps)-(k-1)*ds)/k
  return(c(eps[2]+dh+ds/2,rep(dh+ds,k-2),eps[1]+dh+ds/2)/par("din")[2])
}
pdf_copy <- function (
  width=6,
  height=width*knitr::opts_current$get("fig.height")/knitr::opts_current$get("fig.width"),
  plot.id=NULL,
  filename
) {
  if (missing(filename)) {
    file.id <- if (is.null(plot.id)) { "" } else { paste0("_",plot.id) }
    filename <- knitr::fig_path(paste(file.id,".pdf",sep=""))
  }
  cat("pdf version at:",filename)
  dev.print( file=filename, device=pdf,
             width=width, height=height,
             pointsize=10,
             family="sans")
}
```

Render this, for instance, like:
  ```
templater::render_template("summarize_run.Rmd",output="lostruct_results/type_snp_size_10000_jobid_324902/run_summary.html",change.rootdir=TRUE)
```

```{r data_setup, include=FALSE}
if (!file.exists("config.json")) {
  stop(paste("File", file.path(getwd(), "config.json"), "does not exist. Cannot continue."))
}
opt <- fromJSON("config.json")
if (is.null(opt$weights)) { opt$weights <- 1 }

# original data files
chroms <- opt$chrom_names
bcf.files <- opt$bcf_files
names(bcf.files) <- chroms

sample.ids <- vcf_samples(bcf.files[1])
warning(opt$sample_info)
if (!is.null(opt$sample_info)) {
  samp.file <- opt$sample_info
  samps <- read.table(samp.file,sep="\t",header=TRUE, stringsAsFactors=TRUE)
  names(samps) <- tolower(names(samps))
  # hack for msprime output
  if (all(grepl("^msp_",sample.ids)) & is.numeric(samps$id)) {
    samps$id <- factor(paste0("msp_", samps$id))
  }
  drop_ids <- setdiff(sample.ids, levels(samps$id))
  if (length(drop_ids)>0) {
    warning(sprintf("These samples have no information in the samples file, %s:\n  %s.", samp.file, paste(drop_ids,collapse=" ")))
  }
  samps <- droplevels( samps[match(sample.ids,samps$id),] )
  samps$population <- factor(samps$population)
} else {
  #warning("No population information in the sample file, %s.", samp.file)
  samps <- data.frame(
    ID=sample.ids,
    population=factor(rep("pop",length(sample.ids))) )
}


