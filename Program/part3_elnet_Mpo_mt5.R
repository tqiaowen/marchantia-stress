# GRN Part 3 
# Started: 17 February 2022
# Aims:
#       1) Generate the linear regression models for gene-TF and TF-TF
#       2) Regularisations: LASSO/ elastic net
#       3) Intermediate outputs per response variable
#       4) Final output as a GRN with [goi_name, TF_name, coefficients(beta), gene_name (of TFs), relvar, s, lambda, cv]

# # Working paths (windows)
# wdir <- "G:/My Drive/Projects/Marchantia_2019/grn/3_GRNmodels"
# #odir <- file.path(wdir, "output/elnet")
# mat_path <- file.path(wdir, "all_stress_mt5_nodum.tsv")
# tf_path <- "G:/My Drive/Projects/Marchantia_2019/grn/TF_collate/PlantTFDB_Mpov5r1_prediction_plusTFDB.txt"

# Working paths (ChuckNorris)
ddir <- "/home/qiaowen/Marchantia_2019/dep"
wdir <- "/mnt/md2/qiaowen/Mpo_GRN_models"
odir <- file.path(wdir, "output")
mat_path <- file.path(ddir, "all_stress_mt5_nodum.tsv")
tf_path <- file.path(ddir, "PlantTFDB_Mpov5r1_prediction_plusTFDB.txt")

## Script to create the elastic net derived grn
source('deps/wrap_elnetv3_Mpo.r')
#library(edgeR)
library(comprehenr)
set.seed(2019)
options(stringsAsFactors = FALSE)
# note to self: cd to /home/qiaowen/Marchantia_2019/scripts/GRN_code/
# models calculated so far: [L, lasso]
#mat_type <- "L"
cmode_name <- "lasso"

extract_matrix <- function(mat, x) {
  x_names <- to_vec(for(z in colnames(mat)) if(grepl(x, z, fixed = TRUE) & !grepl('control', z, fixed = TRUE)) z)
  new_matrix <- subset(mat, select = x_names)
  return(new_matrix)
}

#stresses <- c("L", "D", "H", "C", "S","M", "N", "all")
stresses <- c("all")

for (mat_type in stresses){
  print(paste("Building models for subset: ", mat_type))
  # load gene expression matrix and normalise (log transform followed by z transform)
  mat <- as.matrix(read.table(mat_path, header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))
  resdir <- file.path(wdir, cmode_name, mat_type)
  
  if (mat_type != "all"){
    mat <- extract_matrix(mat, mat_type)
  }
  
  # kick out genes that are completely '0' across all conditions (aftifact of subsetting)
  mat <- mat[rowSums(mat) != 0,]
  
  ###
  # Checks if matrix contains zeros
  if (sum(mat>0) > 0){
    print("Warning: '0' present in gene expression matrix.")
    minval <- min(mat[mat > 0])
    print(paste("Minimum expression value:", minval))
    print("Replacing zeros with 1e-12")
    mat[mat == 0] <- 0.000000000001
  }
  
  # Log transforms matrix
  log_mat <- log(mat)
  #sum(colSums(log_mat == -Inf))
  gene_dat <- t(scale(t(log_mat)))
  
  #import TF annotation - names should be the same as rownames of the read data
  TF <- read.delim(tf_path, header = FALSE)
  tfs <- TF[, 1]
  tfs <- tfs[tfs %in% rownames(gene_dat)] # grab only TFs that are in current dataset
  #print(length(tfs))
  
  #Do elastic net regression analysis for each TF using all other TFs as predictors 
  #THIS STEP TAKES VERY LONG - RUN ON SERVER OR wITH MORE THEN 8 THREADS
  elnet_res <- wrap_elnet(gene_dat, resdir=resdir, thrsh=0, tfs=tfs, parallel=64, cmode=cmode_name)
  elnet_all <- elnet_res[,c('Gene.ID', 'predicted', 'rel.coeff')]
  colnames(elnet_all) <- c('from', 'to', 'weight' )
  elnet_all <- elnet_all[order(abs(elnet_all$weight), decreasing = TRUE),]
  save(elnet_all, file=file.path(resdir,'elnet_all.obj'))
  write.table(elnet_all, file.path(resdir, 'elnet_all.txt'), sep = "\t", col.names = NA)
}