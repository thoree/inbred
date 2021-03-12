#'
#' Compare CI-s (nonparametric, parametric bootstrap, and asymptotic)
#'
#' The purpose is to compare confidence intervals based on nonparametric bootstrap,
#' parametric bootstrap, and asymptotic. This is done for the simplest case: iid SNPs
#' with frequency 0.5. We consider pairwise relationships and assume the individuals
#' do not share two alleles IBD. Percentile t, bca
#'
#' @param theta Double kappa0.
#' @param x ped object.
#' @param ids Id of pair.
#' @param n Integer. No of markers.
#' @param N Integer. No of simulations.
#' @param B Integer. No of bootstraps.
#' @param seed Integer.


#' @return Returns table.
#' @export
#' @examples
#' library(forrel)
#' library(pedprobr)
#' x = halfSibPed(1)
#' ids = c(4,5)
#' # ids = c(1,3)
#' n = 100 #markers
#' N = 2 #simulations
#' B = 100 # bootstraps
#' theta = 0.5
#' # theta = 1
#' seed = 17
#' # Frequencies for bootstrap
#' freq = list()
#' for (i in 1:n){
#' freq[[i]] =  list(alleles = 1:2, afreq = c(0.5, 0.5))
#' }
#' x = setMarkers(x, locusAttributes = freq)
#' foo = compareCI(theta, x, ids,freqList, n, N, B, seed)

compareCI <- function(theta, x, ids, n = 2, N = 2, B = 2, seed = NULL){
  # Find genotype probabilities
  set.seed(seed)
  p = findp(x, ids)
  # Simulate from multionomial and find sufficient statistics
  tab = rmultinom(N, n, p)
  rownames(tab) = c("n11", "n12", "n13", "n21", "n22", "n23",
                    "n31", "n32", "n33" )
  n1 = tab["n11",] + tab["n33",]
  n2 = tab["n13",] + tab["n31",]
  n3 = tab["n22",]
  n4 = tab["n12",] + tab["n21",]+tab["n23",] + tab["n32",]

  ## MLE estimate and asymptotic confidence interval
  theta.MLE = pmin(c(2/(1+n1/n2)), 1)
  se = sqrt(1/(n*I1(theta.MLE)))
  CI = matrix(c("lower" = theta.MLE - 1.96*se,
                "upper" = theta.MLE + 1.96*se), ncol = 2)
  inside = CI[,1] <= theta & CI[,2] >= theta
  coverage.MLE = sum(inside)/N
  # Add point estimate and indicator for being in interval
  MLE = data.frame(lower = pmax(CI[,1],0),
                   upper = pmin(CI[,2], 1), theta.MLE, inside)
  MLE.average = c(apply(MLE[,-4],2,mean), coverage =sum(inside)/N)

  # Parametric bootstrap
  CI.para = matrix(ncol = 2, nrow = N)
  theta.para = rep(NA, N)
  for (j in 1:N){
    boot1 = ibdBootstrap(x,kappa = c(theta, 1-theta, 0),  N = B,
                         plot = F)
    CI.para[j,] = quantile(boot1[,1], probs = c(0.025, 0.975))
    theta.para[j] = mean(boot1[,1])
  }
  inside = CI.para[,1] <= theta & round(CI.para[,2],6) >= theta
  coverage.para = sum(inside)/N
  parametric = data.frame(lower = pmax(CI.para[,1],0),
                          upper = pmin(CI.para[,2], 1), theta.para, inside)
  parametric.average = c(apply(parametric[,-4],2,mean), coverage =sum(inside)/N)


  # Nonparametric bootstrap # forrel::nonParamBoot
  CI.non = matrix(ncol = 2, nrow = N)
  theta.non = rep(NA, N)
  CI.summary.non = NA
  for (j in 1:N){
    x1 = profileSim(x,  1, ids, verbose = F)[[1]]
    boot1 = ibdBootstrap(x1, ids, N = B, param = "kappa",
                         freqList = freqList, plot = F, method ="nonparametric")
    CI.non[j,] = quantile(boot1[,1], probs = c(0.025, 0.975))
    theta.non[j] = mean(boot1[,1])
  }
  inside = CI.non[,1] <= theta & round(CI.non[,2],6) >= theta
  coverage.non = sum(inside)/N
  nonparametric = data.frame(lower = pmax(CI.non[,1],0),
                             upper = pmin(CI.non[,2], 1), theta.non, inside)
  nonparametric.average = c(apply(nonparametric[,-4],2,mean), coverage =sum(inside)/N)

  list(MLE = MLE, MLE.average = MLE.average,
       parametric = parametric, parametric.average = parametric.average,
       nonparametric = nonparametric, nonparametric.average = nonparametric.average)
}


findp = function(x, ids){
  x = halfSibPed(1)
  m = marker(x, alleles = c("a", "b"), afreq = c(0.5, 0.5))
  tab0 = oneMarkerDistribution(x, ids = ids,
                               partialmarker = m, verbose = F)
  as.vector(tab0)
}

I1 = function(theta){
  p1.theta = p1(k0 = theta, k1 = 1-theta, k2 = 0)
  p2.theta = p2(k0 = theta)
  (1/256)*((1-p1.theta)/p1.theta + (1-p2.theta)/p2.theta + 2)
}

p1 = function(p = 0.5, k0 = 0.25, k1 = 0.5, k2 = 0.25)
  k0 * p^4 + k1 *p^3 + k2 * p^2

p2 = function(p = 0.5, k0 = 0.25)
  k0 * p^4

