data.dir <- getwd()

## load datasets
clinical.dat<-readRDS(file.path(data.dir, "data_clean/clinical.rds"))
protein.dat<-readRDS(file.path(data.dir, "data_clean/protein.rds"))
mirna.dat<-readRDS(file.path(data.dir, "data_clean/mirna.rds"))
mrna.dat<-readRDS(file.path(data.dir, "data_clean/mrna.rds"))
mutations.dat<-readRDS(file.path(data.dir, "data_clean/mutations.rds"))
methylation.dat<-readRDS(file.path(data.dir, "data_clean/methylation.rds"))

#Install and load packages
library(mlr3verse)
library(mlr3pipelines) 
library(precrec) ## mlr3 plots 
library(rpart) ## classification trees 
library(xgboost) ## xgboost 
library(glmnet) ## ealstic net 
library(future) #parallel running
library(Cairo) #save figure on hpc

plan(multisession,workers=20)


metrics = c("auc", "acc", "sensitivity", "specificity", "precision", "recall", "fbeta")
metrics = paste("classif", metrics, sep = ".")
metrics = lapply(metrics, msr)

# Create Elastic net leaner
elnet = (
  po("imputesample") %>>%
    po("scale") %>>%
    po("pca") %>>%
    lrn(
      "classif.cv_glmnet",
      predict_type="prob"))
elnet = as_learner(elnet)


#Create XGBoost leaner
xgb = (
  po("imputesample") %>>%
    po("scale")%>>%
    po("pca")%>>%
    lrn(
      "classif.xgboost",
      predict_type='prob',
      objective = "binary:logistic",
      eta = 0.3,
      max_depth = 6,
      min_child_weight = 1,
      subsample = 0.5,
      colsample_bytree = 0.5
    ))
xgb=as_learner(xgb)


#Predict breast cancer with protein abundance
protein.bc.dat = data.frame(clinical.dat, t(protein.dat))
task_protein = as_task_classif(
  x=protein.bc.dat,
  target="pfi",
  id="Breast cancer PFI predicted by protein abundance"
)

task_protein$col_roles$stratum = task_protein$target_names

set.seed(43)
dat.parts.protein = partition(task_protein, 3/4)

#Train for elastic net
elnet$train(task_protein, dat.parts.protein$train)

preds.protein.elnet=elnet$predict(task_protein, dat.parts.protein$test)
preds.protein.elnet$confusion

#Train for xgboost
xgb$train(task_protein, dat.parts.protein$train)

preds.xgb.protein=xgb$predict(task_protein, dat.parts.protein$test)
preds.xgb.protein$confusion



#compare performance by benchmarking
design.protein = benchmark_grid(task_protein, c(elnet, xgb), rsmp("cv",folds=3))
bm.protein = benchmark(design.protein)

bm.protein$aggregate(measures = metrics)

png("figures/Boxplot_p.png")
autoplot(bm.protein, measure = msr("classif.auc"))
dev.off()

png("figures/AUC_p.png")
autoplot(bm.protein, type = "roc")
dev.off()

png("figures/PRC_p.png")
autoplot(bm.protein, type = "prc")
dev.off()


##Predict breast cancer by methylation
methylation.bc.dat = data.frame(clinical.dat, t(methylation.dat))
task_meth = as_task_classif(
  x=methylation.bc.dat,
  target="pfi",
  id="Breast cancer PFI predicted by methylation"
)

task_meth$col_roles$stratum = task_meth$target_names

set.seed(47)
dat.parts.meth = partition(task_meth, 3/4)

#Train for elastic net
elnet$train(task_meth, dat.parts.meth$train)

preds.meth.elnet=elnet$predict(task_meth, dat.parts.meth$test)
preds.meth.elnet$confusion

#Train for xgboost
xgb$train(task_meth, dat.parts.meth$train)

preds.xgb.meth=xgb$predict(task_meth, dat.parts.meth$test)
preds.xgb.meth$confusion



#compare performance by benchmarking
design.meth = benchmark_grid(task_meth, c(elnet, xgb), rsmp("cv",folds=3))
bm.meth = benchmark(design.meth)

bm.meth$aggregate(measures = metrics)

png("figures/Boxplot_meth.png")
autoplot(bm.meth, measure = msr("classif.auc"))
dev.off()

png("figures/AUC_meth.png")
autoplot(bm.meth, type = "roc")
dev.off()

png("figures/PRC_meth.png")
autoplot(bm.meth, type = "prc")
dev.off()



##Predict breast cancer by Micro RNA
mirna.bc.dat = data.frame(clinical.dat, t(mirna.dat))
task_mirna = as_task_classif(
  x=mirna.bc.dat,
  target="pfi",
  id="Breast cancer PFI predicted by Micro RNA"
)

task_mirna$col_roles$stratum = task_mirna$target_names

set.seed(93)
dat.parts.mirna = partition(task_mirna, 3/4)

#Train for elastic net
elnet$train(task_mirna, dat.parts.mirna$train)

preds.mirna.elnet=elnet$predict(task_mirna, dat.parts.mirna$test)
preds.mirna.elnet$confusion

#Train for xgboost
xgb$train(task_mirna, dat.parts.mirna$train)

preds.xgb.mirna=xgb$predict(task_mirna, dat.parts.mirna$test)
preds.xgb.mirna$confusion



#compare performance by benchmarking
design.mirna = benchmark_grid(task_mirna, c(elnet, xgb), rsmp("cv",folds=3))
bm.mirna = benchmark(design.mirna)

bm.mirna$aggregate(measures = metrics)

png("figures/Boxplot_mirna.png")
autoplot(bm.mirna, measure = msr("classif.auc"))
dev.off()

png("figures/AUC_mirna.png")
autoplot(bm.mirna, type = "roc")
dev.off()

png("figures/PRC_mirna.png")
autoplot(bm.mirna, type = "prc")
dev.off()




## Predict breast cancer by Messenger RNA
mrna.bc.dat = data.frame(clinical.dat, t(mrna.dat))
task_mrna = as_task_classif(
  x=mrna.bc.dat,
  target="pfi",
  id="Breast cancer PFI predicted by Messenger RNA"
)

task_mrna$col_roles$stratum = task_mrna$target_names

set.seed(63)
dat.parts.mrna = partition(task_mrna, 3/4)

#Train for elastic net
elnet$train(task_mrna, dat.parts.mrna$train)

preds.mrna.elnet=elnet$predict(task_mrna, dat.parts.mrna$test)
preds.mrna.elnet$confusion

#Train for xgboost
xgb$train(task_mrna, dat.parts.mrna$train)

preds.xgb.mrna=xgb$predict(task_mrna, dat.parts.mrna$test)
preds.xgb.mrna$confusion



#compare performance by benchmarking
design.mrna = benchmark_grid(task_mrna, c(elnet, xgb), rsmp("cv",folds=3))
bm.mrna = benchmark(design.mrna)

bm.mrna$aggregate(measures = metrics)

png("figures/Boxplot_mrna.png")
autoplot(bm.mrna, measure = msr("classif.auc"))
dev.off()

png("figures/AUC_mrna.png")
autoplot(bm.mrna, type = "roc")
dev.off()

png("figures/PRC_mrna.png")
autoplot(bm.mrna, type = "prc")
dev.off()




## Predict breast cancer by mutations
mutations.bc.dat = data.frame(clinical.dat, t(mutations.dat))
task_mu = as_task_classif(
  x=mutations.bc.dat,
  target="pfi",
  id="Breast cancer PFI predicted by mutations"
)

task_mu$col_roles$stratum = task_mu$target_names

set.seed(88)
dat.parts.mu = partition(task_mu, 3/4)

#Train for elastic net
elnet$train(task_mu, dat.parts.mu$train)

preds.mu.elnet=elnet$predict(task_mu, dat.parts.mu$test)
preds.mu.elnet$confusion

#Train for xgboost
xgb$train(task_mu, dat.parts.mu$train)

preds.xgb.mu=xgb$predict(task_mu, dat.parts.mu$test)
preds.xgb.mu$confusion



#compare performance by benchmarking
design.mu = benchmark_grid(task_mu, c(elnet, xgb), rsmp("cv",folds=3))
bm.mu = benchmark(design.mu)

bm.mu$aggregate(measures = metrics)

png("figures/Boxplot_mu.png")
autoplot(bm.mu, measure = msr("classif.auc"))
dev.off()

png("figures/AUC_mu.png")
autoplot(bm.mu, type = "roc")
dev.off()

png("figures/PRC_mu.png")
autoplot(bm.mu, type = "prc")
dev.off()