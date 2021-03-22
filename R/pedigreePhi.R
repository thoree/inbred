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
#' @param plot Logical
#'
#' @return Returns a l
#'
#'
#' @return Returns a list with two elements, a dataframe and the result
#'   of last simulation called `lastBoot`(see documentation of `ibdBootstrap`).
#'   The columns of the data frame are
#'
#'   * `phi.hat` The estimate for each bootstrap
#'
#'   * `skew` See details
#'
#'   * `dist` See documentation of `ibdBootstrap`
#'
#'   * `lower` Lower bound for interval
#'
#'
#'   * `upper` Upper bound for interval
#'
#'   * `coverage`

#' @details See note coverage.pdf
#'
#' @export
#'
#' @examples
#' library(forrel)
#' library(pedprobr)
#' library(ribd)
#' library(moments)
#'
#' # Example Estimate coverage and more
#' n = 100 # no of markers
#' p = rep(0.5, n)
#' freq = list()
#' for (i in 1:n)
#'   freq[[i]] =  list(afreq = c("1" = p[i], "2" = 1- p[i]))
#' ped = quadHalfFirstCousins()
#' ped = setMarkers(ped, locusAttributes = freq)
#' N = 2 # no of confidence intervals
#' B = 100 # no of  bootstraps
#' seed = 2
#' kappa = ribd::kappaIBD(ped, ids = leaves(ped))
#' res1 = pedigreePhi(kappa = kappa, ped = ped, ids,
#'                    N = N, B = B, seed = seed, level = 0.95, plot = T)
#' round(res1$tab[c(1:2,N+1),], 4)
#' boot1 = res1$lastBoot
#' phi.hat = 0.25*boot1$k1 + 0.5*boot1$k2
#' plot(density(phi.hat), main = "", xlab = "Kinship coefficient")
#' qqnorm(phi.hat)


pedigreePhi <- function(kappa, ped, ids = NULL, N = 2, B = 2,
                        seed = NULL, level = 0.95, plot = F){
  if(!is.null(seed))
    set.seed(seed)
  CI = matrix(ncol = 2, nrow = N)
  phi = 0.25*kappa[2] + 0.5*kappa[3]
  phi.hat = skew = dist = rep(NA, N)

  for (j in 1:N){
    boot1 = ibdBootstrap(ped, kappa = kappa,  N = B, plot = plot)
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
  list(tab = as.data.frame(res2), lastBoot = boot1)
}

ci = function(x, level = 0.95){
  est = mean(x)
  n = length(x)
  se = sd(x)/sqrt(n)
  z = qnorm((1-level)/2,lower.tail = FALSE)
  CI = c(est - z*se, est + z*se)
  CI
}



