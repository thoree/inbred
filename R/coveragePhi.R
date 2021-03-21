#'
#' Coverage for bootstrap confidence intervals
#'
#' Currently we only consider parametric bootstrap and the kinship coefficient.
#'
#' @param kappa Probability vector of length 3.
#' @param ped ped object with allele frequencies.
#' @param ids Id of pair.
#' @param N Integer. No of simulations.
#' @param B Integer. No of bootstraps.
#' @param seed Integer.
#' @param level Double
#' @return Returns a l
#'
#' @details The equations
#'
#' @export
#' @examples
#' library(forrel)
#' library(pedprobr)
#' library(ribd)
#' library(moments)
#'
#' # Example 1a Estimate coverage and more
#' n = 1000 # no of markers
#' p = rep(0.5, n)
#' freq = list()
#' for (i in 1:n)
#'   freq[[i]] =  list(afreq = c("1" = p[i], "2" = 1- p[i]))
#' ped = quadHalfFirstCousins()
#' ped = setMarkers(ped, locusAttributes = freq)
#' N = 400 # no of confidence intervals
#' B = 400 # no of  bootstraps
#' seed = 2
#' kappa = ribd::kappaIBD(ped, ids = leaves(ped))
#' kappa = c(0.95, 0.05, 0)
#' res1 = coveragePhi(kappa = kappa, ped = ped, ids, N = N, B = B, seed = seed, level = 0.95)
#' round(res1[c(1:2,N+1),], 4)
#'
#' Example 1b Plotting for Example 1a
#' windows()
#' par(mfrow = c(1,3))
#' boot1 = ibdBootstrap(ped, kappa = kappa,  N = B, plot = T)
#' phi.hat = 0.25*boot1$k1 + 0.5*boot1$k2
#' plot(density(phi.hat), main = "", xlab = "Kinship coefficient")
#' qqnorm(phi.hat)


coveragePhi <- function(kappa, ped, ids = NULL, N = 2, B = 2, seed = NULL, level = 0.95){
  if(!is.null(seed))
    set.seed(seed)
  CI = matrix(ncol = 2, nrow = N)
  phi = 0.25*kappa[2] + 0.5*kappa[3]
  phi.hat = skew = dist = rep(NA, N)

  for (j in 1:N){
    boot1 = ibdBootstrap(ped, kappa = kappa,  N = B, plot = F)
    phis = 0.25*boot1[,2] + 0.5*boot1[,3]
    phi.hat[j] = mean(phis)
    skew[j] = skewness(phis)
    dist[j]= mean(boot1$dist)
    CI[j,] = ci(phis, level = level)
   }
  coverage = (CI[,1] <= phi) & (CI[,2] >= phi)
  res1 = cbind(phi.hat, skew, dist, lower = CI[,1], upper = CI[,2], coverage)
  res2 = rbind(res1, apply(res1, 2, mean))
  rownames(res2) = c(paste0("sim", 1:N), "average")
  as.data.frame(res2)
}

ci = function(x, level = 0.95){
  est = mean(x)
  n = length(x)
  se = sd(x)/sqrt(n)
  z = qnorm((1-level)/2,lower.tail = FALSE)
  CI = c(est - z*se, est + z*se)
  CI
}



