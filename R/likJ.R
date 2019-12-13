#' Likelihood for two Possibly Inbred Individuals as a function of
#' the Condensed 'Jacquard' Coefficients
#'
#' @param a Vector of positive integers (allele 1 of individual 1)
#' @param b Vector of positive integers (allele 2 of individual 1)
#' @param cc Vector of positive integers (allele 1 of individual 2)
#' @param d Vector of positive integers (allele 2 of individual 2)
#' @param pa Double vector of allele frequencies
#' @param pb Double vector of allele frequencies
#' @param pc Double vector of allele frequencies
#' @param pd Double vector of allele frequencies
#' @param Delta Double vector of length 9 summing to unity
#' @export
#' @examples
#' p = c("1" = 0.1, "2" = 0.9)
#' a  = c(rep(1,6), rep(2,3))
#' b  = c(rep(1,3), rep(2,6))
#' cc = rep(c(1,1,2), 3)
#' d  = rep(c(1,2,2), 3)
#' pa = p[a]
#' pb = p[b]
#' pc = p[cc]
#' pd = p[d]
#' Delta = c(0, 0, 0, 0, 0, 0, 1,0, 0)
#' l2 = likJ(a,b,cc,d, pa,pb,pc,pd, Delta = Delta)
#' sum(l2) == 1
likJ = function(a,b,cc,d, pa, pb, pc, pd, Delta){
  I12 = a==b
  I13 = a==cc
  I14 = a==d
  I23 = b==cc
  I24 = b==d
  I34 = cc==d
  nMarkers = length(a)

  # create objects of correct size
  PG = matrix(ncol=9, nrow=nMarkers)
  case = character(nMarkers)

  index1 = I12 & I34 & I13 # case 1
  index2  = I12 & I34 & !I13 # case 2
  index3a = I12 & I13 & !I14 # case 3a
  index3b = I12 & !I13 & I14 # case 3b
  index4 = I12 & !I13 & !I14 & !I34 # case = 4
  index5a = I34 & I13 & !I12 # case 5a
  index5b = I34 & I23 & !I12  # case 5b
  index6 = I34 & !I12 & !I13 & !I23 # case 6
  index7a = I13 & I24 & !I12  # case 7a
  index7b = I14 & I23 & !I12  # case 7b
  index8a = I13 & !I24 & !I12 & !I34 # case 8a
  index8b = I14 & !I23 & !I12 & !I34 # case 8b
  index8c = I23 & !I14 & !I12 & !I34 # case 8c
  index8d = I24 & !I13 & !I12 & !I34 # case 8d
  index9 = !I12 & !I34 & !I13 & !I24 & !I14 & !I23 # case 9

  if(any(index1)) {
    case[index1] = '1'
    p = pa[index1]
    PG[index1,] = cbind(p,p^2,p^2,p^3,p^2,p^3,p^2, p^3,p^4)
  }

  if(any(index2)) {
    case[index2] = '2'
    p = pa[index2]
    q = pc[index2]
    PG[index2,] = cbind(0,p*q,0,p*q^2,0,p^2*q,0,0,p^2*q^2)
  }

  if(any(index3a)) {
    case[index3a] = '3a'
    p = pa[index3a]
    q = pd[index3a]
    PG[index3a,] = cbind(0,0,p*q,2*p^2*q,0,0,0,p^2*q,2*p^3*q)
  }

  if(any(index3b)) {
    case[index3b] = '3b'
    p = pa[index3b]
    q = pc[index3b]
    PG[index3b,] = cbind(0,0,p*q,2*p^2*q,0,0,0,p^2*q,2*p^3*q)
  }

  if(any(index4)) {
    case[index4] = '4'
    p = pa[index4]
    q = pc[index4]
    r = pd[index4]
    PG[index4,] = cbind(0,0,0,2*p*q*r,0,0,0,0,2*p^2*q*r)
  }

  if(any(index5a)) {
    case[index5a] = '5a'
    p  = pc[index5a]
    q =  pb[index5a]
    PG[index5a,] = cbind(0,0,0,0,p*q,2*p^2*q,0,p^2*q,2*p^3*q)
  }

  if(any(index5b)) {
    case[index5b] = '5b'
    p  = pc[index5b]
    q =  pa[index5b]
    PG[index5b,] = cbind(0,0,0,0,p*q,2*p^2*q,0,p^2*q,2*p^3*q)
  }

  if(any(index6)) {
    case[index6] = '6'
    p = pc[index6]
    q = pa[index6]
    r = pb[index6]
    PG[index6,] = cbind(0,0,0,0,0,2*p*q*r,0,0,2*p^2*q*r)
  }

  if(any(index7a)) {
    case[index7a] = '7a'
    p = pa[index7a]
    q = pb[index7a]
    PG[index7a,] = cbind(0,0,0,0,0,0,2*p*q,p*q*(p+q),4*p^2*q^2)
  }

  if(any(index7b)) {
    case[index7b] = '7b'
    p = pa[index7b]
    q = pb[index7b]
    PG[index7b,] = cbind(0,0,0,0,0,0,2*p*q,p*q*(p+q),4*p^2*q^2)
  }

  if(any(index8a)) {
    case[index8a] = '8a'
    p = pa[index8a]
    q = pb[index8a]
    r = pd[index8a]
    PG[index8a,] = cbind(0,0,0,0,0,0,0,p*q*r,4*p^2*q*r)
  }

  if(any(index8b)) {
    case[index8b] = '8b'
    p = pa[index8b]
    q = pb[index8b]
    r = pc[index8b]
    PG[index8b,] = cbind(0,0,0,0,0,0,0,p*q*r,4*p^2*q*r)
  }

  if(any(index8c)) {
    case[index8c] = '8c'
    p = pa[index8c]
    q = pb[index8c]
    r = pd[index8c]
    PG[index8c,] = cbind(0,0,0,0,0,0,0,p*q*r,4*q^2*p*r)
  }

  if(any(index8d)) {
    case[index8d] = '8d'
    p = pa[index8d]
    q = pb[index8d]
    r = pc[index8d]
    PG[index8d,] = cbind(0,0,0,0,0,0,0,p*q*r,4*q^2*p*r)
  }

  if(any(index9)) {
    case[index9] = '9'
    p = pa[index9]
    q = pb[index9]
    r = pc[index9]
    s = pd[index9]
    PG[index9,] = cbind(0,0,0,0,0,0,0,0,4*p*q*r*s)
  }

  lik = as.numeric(PG %*% Delta)
  lik
}


