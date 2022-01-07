#!/usr/bin/Rscript
# Copyright (c) 2013, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
# BEGINCONF
# [pvalue]
# comment=p-value from MetaSKAT_wZ method
# [sample.size]
# # Adjust this for your data !!
# n = 3
# # Adjust this for your data !!
# name = n.pop1, n.pop2, n.pop3
# comment=sample size per group
# ENDCONF
suppressMessages(library(MetaSKAT))
MetaSKAT.VAT <- function (dat, out_type, group_colname,
                          r.corr = 0,
                          pval.method = "optimal",
                          combined.weight = TRUE,
                          is.separate = FALSE,
                          weights.beta = c(1,25),
                          Group_Idx = NULL) {
  # see how many variants we have in the problem
  if (ncol(dat@X) == 1) {
    write("Only one variant found, not a rare variant analysis problem", stderr())
    q("no", 1, FALSE)
  }
  # n.g a numeric value of the number of cohort
  groups <- dat@Y[, group_colname]
  ugroups <- unique(groups)
  n.g <- length(ugroups)
  if (n.g == 1) {
    write("Only one cohort found, not a meta analysis problem", stderr())
    q("no", 1, FALSE)
  }
  # Group_Idx is a vector of group index.
  # Suppose the ﬁrst two cohorts are European cohorts and the
  # last cohort is an African American cohort.
  # If you want to run MetaSKAT with assuming ancestry
  # group speciﬁc heterogeneity, you can set Group_Idx=c(1,1,2),
  # which indicates the ﬁrst two cohorts belong to the same group.
  if (is.null(Group_Idx)) Group_Idx <- seq(1:n.g)
  # Drop the group label
  Y <- subset(dat@Y, select = -c(get(group_colname)))
  m <- ncol(Y)
  # y.list is a list object of phenotypes.
  # It has m elements for m studies.
  # Each element is a vector of phenotypes.
  y.list <- NULL; for (i in 1:n.g) y.list[[i]] <- Y[which(groups == ugroups[i]), m]
  # x.list is a list object of covariates matrices.
  # It has m elements for m studies.
  # Each element is a matrix of covariates
  x.list <- NULL; for (i in 1:n.g) x.list[[i]] <- as.matrix(Y[which(groups == ugroups[i]), -m])
  # Z is matrix of the gene
  # The rows of Z should be matched with phenotypes and covariates.
  # If there are 3 studies, and study 1,2, and 3 have n1, n2, and n3 samples,
  # the ﬁrst n1, n2, and n3 rows of Z should be the genotypes of
  # the ﬁrst, second, and third studies, respectively.
  Z <- NULL; for (i in 1:n.g) Z <- rbind(Z, as.matrix(dat@X[which(groups == ugroups[i]), ]))
  #
  obj <- Meta_Null_Model(y.list, x.list, n.cohort = n.g, out_type = out_type)
  res <- MetaSKAT_wZ(Z, obj, r.corr = r.corr, method = pval.method,
                     combined.weight = combined.weight, weights.beta = weights.beta,
                     is.separate = is.separate,
                     Group_Idx = Group_Idx)
  res.formatted <- list(sample.size = getGroupSampleSize(groups, ugroups), pvalue = res$p.value)
  return(res.formatted)
}

getGroupSampleSize <- function(groups, ugroups) {
  size <- vector()
  for (i in 1:length(ugroups)) {
    size[i] <- length(groups[which(groups == ugroups[i])])
  }
  return(size)
}
