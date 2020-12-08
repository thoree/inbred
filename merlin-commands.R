library(pedprobr)
#> Loading required package: pedtools

# Pedigree and markers
x = halfCousinPed(1)
x = setMarkers(x, list(
  marker(x, "5" = 1, "9" = 1, alleles = 1:2, afreq = c(0.2, 0.8), name = "L1"),
  marker(x, "5" = 1, "9" = 1, alleles = 1:2, afreq = c(0.2, 0.8), name = "L2")))

res = merlin(x, options = '--simul --reruns 2 --save -r 123')
cat(res, sep = "\n")

# Likelihood computation works as before:
likelihoodMerlin(x, markers = 1, verbose = F)
likelihood(x, 1) # same except merlin's rounding



library(pedprobr)
x = swapSex(halfCousinPed(1), c(1, 3, 7))
x = setMarkers(x, list(
  marker(x, "5" = 1, "9" = 1, alleles = 1:2, afreq = c(0.2, 0.8),
         chrom = "X", name = "L1"),
  marker(x, "5" = 1, "9" = 1, alleles = 1:2, afreq = c(0.2, 0.8),
         chrom = "X", name = "L2")))
likelihoodMerlin(x, markers = 1,verbose = F)
likelihood(x,1)
# Formula
p = 0.2
k = ribd::kappaIBD(x,leaves(x))
k[1]*p^2+k[2]*p
