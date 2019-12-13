#' Case determination for a pair of individuals
#'
#' This is a reduced version of the function `likJ`. It only returns
#' the 'casew' variable oflikJ
#'
#' @param a Vector of positive integers (allele 1 of individual 1)
#' @param b Vector of positive integers (allele 2 of individual 1)
#' @param cc Vector of positive integers (allele 1 of individual 2)
#' @param d Vector of positive integers (allele 2 of individual 2)

#' @export
#' @return `case`
#' @examples
#' # Below I go through some cases trying to figure out
#' # what the LD probabilities should be
#'
#' p = c(1,2,3,4,5)/sum(1:5)
#' H = p%o%p # Haplotypeprobs for LE
#'
#' ### Case  1,1
#' a  = c(1, 2)
#' b  = c(1, 2)
#' cc = c(1, 2)
#' d  = c(1, 2)
#' (Case = case(a,b,cc,d))
#' pa = p[a]
#' pb = p[b]
#' pc = p[cc]
#' pd = p[d]
#' Delta = c(1, 0, 0, 0, 0, 0, 0, 0, 0)
#' lik.LE = prod(likJ(a, b, cc, d, pa, pb, pc, pd, Delta))
#' lik.LD = H[a[1], a[2]]
#' lik.LD = geno1(cc[1], d[1], cc[2], d[2], H)
#' abs(lik.LE -lik.LD) < 1e-15
#'
#' ### Case 1,2
#' a  = c(1, 1)
#' b  = c(1, 1)
#' cc = c(1, 2)
#' d  = c(1, 2)
#' (Case = case(a,b,cc,d))
#' pa = p[a]
#' pb = p[b]
#' pc = p[cc]
#' pd = p[d]
#' Delta1 = c(1, 0, 0, 0, 0, 0, 0, 0, 0)
#' Delta2 = c(0, 1, 0, 0, 0, 0, 0, 0, 0)
#' lik.LE = likJ(a[1], b[1], cc[1], d[1], pa[1], pb[1], pc[1], pd[1], Delta1)*
#'          likJ(a[2], b[2], cc[2], d[2], pa[2], pb[2], pc[2], pd[2], Delta2)
#' # P(G1,G2|Jst)= P(G1|Jst, G2)  * P(G2|Jst) = p[a[[1]]]*H[cc[1],cc[2]]
#' lik.LD = p[a[[2]]]*H[cc[1],cc[2]]
#' lik.LD = p[a[[2]]]*geno1(cc[1], d[1], cc[2], d[2], H)
#' abs(lik.LE -lik.LD) < 1e-15
#'
#' ### Case 2,1

#' a  = c(1, 1)
#' b  = c(1, 1)
#' cc = c(2, 1)
#' d  = c(2, 1)
#' (Case = case(a,b,cc,d))
#' pa = p[a]
#' pb = p[b]
#' pc = p[cc]
#' pd = p[d]
#' Delta1 = c(0, 1, 0, 0, 0, 0, 0, 0, 0)
#' Delta2 = c(1, 0, 0, 0, 0, 0, 0, 0, 0)
#' lik.LE = likJ(a[1], b[1], cc[1], d[1], pa[1], pb[1], pc[1], pd[1], Delta1)*
#'          likJ(a[2], b[2], cc[2], d[2], pa[2], pb[2], pc[2], pd[2], Delta2)
#' # P(G1,G2|Jst)= P(G1|Jst, G2)  * P(G2|Jst) = p[a[[1]]]*H[cc[1],cc[2]]
#' lik.LD = p[a[[1]]]*H[cc[1],cc[2]]
#' lik.LD = p[a[[1]]]*geno1(cc[1], d[1], cc[2], d[2], H)
#' abs(lik.LE -lik.LD) < 1e-15
#'
case = function(a,b,cc,d){
  I12 = a==b
  I13 = a==cc
  I14 = a==d
  I23 = b==cc
  I24 = b==d
  I34 = cc==d
  nMarkers = length(a)

  # create objects of correct size
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

  if(any(index1))
    case[index1] = '1'
  if(any(index2))
    case[index2] = '2'
  if(any(index3a))
    case[index3a] = '3a'
  if(any(index3b))
    case[index3b] = '3b'
  if(any(index4))
    case[index4] = '4'
  if(any(index5a))
    case[index5a] = '5a'
  if(any(index5b))
    case[index5b] = '5b'
  if(any(index6))
    case[index6] = '6'
  if(any(index7a))
    case[index7a] = '7a'
  if(any(index7b))
    case[index7b] = '7b'
  if(any(index8a))
    case[index8a] = '8a'
  if(any(index8b))
    case[index8b] = '8b'
  if(any(index8c))
    case[index8c] = '8c'
  if(any(index8d))
    case[index8d] = '8d'
  if(any(index9))
    case[index9] = '9'
  case
}


