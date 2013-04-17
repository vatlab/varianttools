#!/usr/bin/Rscript
# Copyright (c) 2013, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
# Meta SKAT Phenotype From File (MetaSKAT.PFF)
# BEGINCONF
# [pvalue]
# comment=p-value from MetaSKAT_wZ method
# [sample.size]
# # Adjust this for your data !!
# n = 2 
# # Adjust this for your data !!
# name = n.pop1, n.pop2 
# comment=sample size per group
# ENDCONF
# Usage:
# vtools associate variant sample_id \
# -m 'RTest /path/to/this/script.R --name MetaSKAT --phenotype_files "c("file1.txt","file2.txt")" --phenotype_colname "QT_TRAIT" --out_type "C" ' \
# --group_by name2 -j8 --to_db MetaSKAT > result.txt 
suppressMessages(library(MetaSKAT))
MetaSKAT.PFF.VAT <- function (dat, out_type,
                          phenotype_colname,
                          phenotype_files,
                          sample_colname = "sample_name",
                          r.corr = 0,
                          pval.method = "optimal",
                          combined.weight = TRUE,
                          is.separate = FALSE) {
  if (ncol(dat@X) == 1) {
    write("Only one variant found, not a rare variant analysis problem", stderr())
    q("no", 1, FALSE)
  }
  n.g <- length(phenotype_files)
  if (n.g == 1) {
    write("Only one cohort found, not a meta analysis problem", stderr())
    q("no", 1, FALSE)
  }
  #
  # step 1: load data
  #
  phenotype <- NULL
  for (i in 1:n.g) phenotype[[i]] <- read.table(phenotype_files[i], header = T)
  #
  # step 2: determine group label
  #
  Group_Idx <- seq(1:n.g)
  #
  # step 3: match sample names for each group with input genotype data
  #
  valid_rows <- NULL
  sample_size <- vector()
  for (i in 1:n.g) {
    valid_rows[[i]] <- phenotype[[i]][, sample_colname] %in% rownames(dat@X)
    sample_size[i] <- sum(valid_rows[[i]])
  }
  #
  # step 4: create y.list, x.list (covariates: phenotype name and sample name columns are dropped) and Z
  #
  y.list <- NULL; for (i in 1:n.g) {
    y.list[[i]] <- phenotype[[i]][valid_rows[[i]], phenotype_colname]
    if (is.null(y.list[[i]])) {
      write(paste("No valid sample found for group", i), stderr())
      q("no", 1, FALSE)
    }
    if (sum(is.na(y.list[[i]])) < length(y.list[[i]])) {
      write("Missing data in phenotype is not allowed", stderr())
      q("no", 1, FALSE)
    }
  }
  x.list <- NULL; for (i in 1:n.g) {
    x.list[[i]] <- as.matrix(
      phenotype[[i]][valid_rows[[i]], -which(names(phenotype[[i]]) %in% c(phenotype_colname, sample_colname))]
      )
  }
  Z <- NULL; for (i in 1:n.g) {
    Z <- rbind(Z, as.matrix(dat@X[which(rownames(dat@X) %in% phenotype[[i]][, sample_colname]), ]))
  }
  #
  # step 5: run meta SKAT
  #
  obj <- Meta_Null_Model(y.list, x.list, n.cohort = n.g, out_type = out_type)
  res <- MetaSKAT_wZ(Z, obj, r.corr = r.corr, method = pval.method,
                     combined.weight = combined.weight, is.separate = is.separate,
                     Group_Idx = Group_Idx)
  return(list(sample.size = sample_size,
              pvalue = res$p.value))
}

getGroupSampleSize <- function(sample_size) {
  size <- vector()
  for (i in 1:length(sample_size)) {
    size[i] <- paste('GP.', i, '=', sample_size[i], sep = '')
  }
  return(paste(size, collapse=";"))
}
