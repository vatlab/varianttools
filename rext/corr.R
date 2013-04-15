#
# A toy example of the R extension:
# calculate correlation or covariance of phenotype Y and counts of mutations per persion
#
corr = function (dat, statistic = "cor") {
  y = dat@Y[,ncol(dat@Y)]
  x = apply(dat@X, 1, sum)
  return (list(sample.size=c(length(y), "integer", "sample size"),
               stat=c(get(statistic)(x, y), 'float', 'statistic')))
}
