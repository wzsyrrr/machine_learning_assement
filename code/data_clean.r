data.dir <- getwd()

##Downloading data from google drive


library(data.table) ## for loading the data 

## Function for loading tab-delimited spreadsheets
read.dataset=function(filename, ...) {
  cat("reading", basename(filename), "... ")
  ## read in tab-delimited spreadsheet
  x=fread(
    filename,
    header=T,
    stringsAsFactors=F,
    sep="\t",
    check.names=F,
    ...)
  ## remove any duplicate rows (identified by the first column)
  x=x[match(unique(x[[1]]), x[[1]]),]
  ## make the first column the rownames of the data frame
  x=data.frame(x,row.names=1,stringsAsFactors=F,check.names=F)
  cat(nrow(x), "x", ncol(x), "\n")
  x
}

## The format of sample identifiers/barcodes is described here:
## https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
##
## Function extracts the participant identifier from a sample id/barcode.
extract.participant <- function(id)
  sub("TCGA-[^-]+-([^-]+)-.*", "\\1", id)

## load datasets
clinical.dat=read.dataset(file.path(data.dir, "data_raw/clinical.txt"))
protein.dat=read.dataset(file.path(data.dir, "data_raw/protein.txt"))
mirna.dat=read.dataset(file.path(data.dir, "data_raw/mirna.txt"))
mrna.dat=read.dataset(file.path(data.dir, "data_raw/mrna.txt"))
mutations.dat=read.dataset(file.path(data.dir, "data_raw/mutations.txt"))
methylation.dat=read.dataset(file.path(data.dir, "data_raw/methylation.txt"))


## harmonize datasets
protein.ids=extract.participant(colnames(protein.dat))
clinical.ids=rownames(clinical.dat)
mirna.ids=colnames(mirna.dat)
mutations.ids=colnames(mutations.dat)
methylation.ids=extract.participant(colnames(methylation.dat))
mrna.ids=extract.participant(colnames(mrna.dat))

common.ids=intersect(clinical.ids, protein.ids)
common.ids=intersect(common.ids, mirna.ids)
common.ids=intersect(common.ids, mutations.ids)
common.ids=intersect(common.ids, methylation.ids)
common.ids=intersect(common.ids, mrna.ids)

clinical.dat=clinical.dat[match(common.ids, clinical.ids),]
protein.dat=protein.dat[,match(common.ids, protein.ids)]
mirna.dat=mirna.dat[,match(common.ids, mirna.ids)]
mrna.dat=mrna.dat[,match(common.ids, mrna.ids)]
methylation.dat=methylation.dat[,match(common.ids, methylation.ids)]
mutations.dat=mutations.dat[,match(common.ids, mutations.ids)]

##Find all columns with constatnt value and remove them
constant_columns <- sapply(mirna.dat, function(x) length(unique(x)) == 1)
mirna.dat <- mirna.dat[, !constant_columns]


## restrict clinical dataset to clinical variables of interest
clinical.vars=c(
  "age.at.diagnosis","estrogen.receptor.status",
  "progesterone.receptor.status",
  "lymphocyte.infiltration","necrosis.percent")
target.var="pfi" ## outcome variable

clinical.dat=clinical.dat[,c(target.var,clinical.vars)]
clinical.dat$estrogen.receptor.status=ifelse(clinical.dat$estrogen.receptor.status=="positive",1,0)
clinical.dat$progesterone.receptor.status=ifelse(clinical.dat$progesterone.receptor.status=="positive",1,0)


convert_column <- function(x) {
  # Remove NA for checking unique values
  unique_values <- unique(na.omit(x))
  # Check if column contains only 1s and 0s
  if (all(unique_values %in% c(0, 1)) && is.numeric(x)) {
    return(as.logical(x))
  } else {
    return(as.numeric(x))
  }
}

clinical.dat[] <- lapply(clinical.dat, convert_column)


## remove features with > 20% missing values
meth.missing.pct=sapply(methylation.dat, function(v) mean(is.na(v)))
methylation.dat=methylation.dat[,meth.missing.pct < 0.2]
mirna.missing.pct=sapply(mirna.dat, function(v) mean(is.na(v)))
mirna.dat=mirna.dat[,mirna.missing.pct < 0.2]
mrna.missing.pct=sapply(mrna.dat, function(v) mean(is.na(v)))
mrna.dat=mrna.dat[,mrna.missing.pct < 0.2]
mu.missing.pct=sapply(mutations.dat, function(v) mean(is.na(v)))
mutations.dat=mutations.dat[,mu.missing.pct < 0.2]
protein.missing.pct=sapply(protein.dat, function(v) mean(is.na(v)))
protein.dat=protein.dat[,protein.missing.pct < 0.2]


writeRDS(methylation.dat,file=file.path(data.dir, "data_clean/methylation.rds"))
writeRDS(mrna.dat,file=file.path(data.dir, "data_clean/methylation.rds"))
writeRDS(mirna.dat,file=file.path(data.dir, "data_clean/methylation.rds"))
writeRDS(mutations.dat,file=file.path(data.dir, "data_clean/methylation.rds"))
writeRDS(protein.dat,file=file.path(data.dir, "data_clean/methylation.rds"))
writeRDS(clinical.dat,file=file.path(data.dir, "data_clean/methylation.rds"))