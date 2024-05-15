data.dir <- getwd()

## load datasets
clinical.dat<-readRDS(file.path(data.dir, "data_clean/clinical.rds"))
protein.dat<-readRDS(file.path(data.dir, "data_clean/protein.rds"))
mirna.dat<-readRDS(file.path(data.dir, "data_clean/mirna.rds"))
mrna.dat<-readRDS(file.path(data.dir, "data_clean/mrna.rds"))
mutations.dat<-readRDS(file.path(data.dir, "data_clean/mutations.rds"))
methylation.dat<-readRDS(file.path(data.dir, "data_clean/methylation.rds"))


omics = list(
  protein.dat=protein.dat,
  mirna.dat=mirna.dat,
  mrna.dat=mrna.dat,
  mutations.dat=mutations.dat,
  methylation.dat=methylation.dat)


#Transpose omics matrices
for (i in 1:length(omics)) {
  omics[[i]]=t(omics[[i]])
}

#Centre and scale
for (i in 1:length(omics)) {
  omics[[i]]=scale(omics[[i]], center = T, scale = T)
}

#Remove columns with more than 20% miising values
for (i in 1:length(omics)) {
  # Calculate the percentage of missing values for each column
  pct_missing <- apply(omics[[i]], 2, function(x) mean(is.na(x)))
  
  # Identify columns where more than 20% of values are missing
  cols_to_remove <- pct_missing > 0.20
  
  # Report the number of columns being removed
  if (any(cols_to_remove)) {
    cat(sum(cols_to_remove), "columns with more than 20% missing values are being removed from", names(omics)[i], "\n")
    
    # Remove the columns with more than 20% missing values
    omics[[i]] <- omics[[i]][, !cols_to_remove]
  }
}

library(mlr3verse)
library(mlr3pipelines) 
library(bcv)


omics <- lapply(omics, function(df) {
  # Check if there are any missing values in the dataframe
  if (any(is.na(df))) {
    # Apply impute.svd
    imputation_result <- impute.svd(df)
    # Return the imputed data
    return(imputation_result$x)
  } else {
    # Return the original data if no imputation is needed
    return(df)
  }
})
saveRDS(imputed_omics_list, file = file.path(data.dir,"imp_dat.rds"))

min.prop = 0.8
r = sapply(omics, function(omics) {
  sdev = prcomp(omics)$sdev ## std dev of each principal component
  which(cumsum(sdev^2/sum(sdev^2)) > min.prop)[1] })

r
cbind(
  features=sapply(omics, ncol),
  "rank estimate"=r)
