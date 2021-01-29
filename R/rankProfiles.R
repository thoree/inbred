#' Ranking profiles of pedigree member
#'
#' @param A `ped` object with attached markers
#' @param id Name of individual to be predicted
#' @param Names Namesor indices of the markers to be included. By default, all markers.
#' @param maxPerMarker No of top candidates per marker, default all
#' @param A logical, by default TRUE
#' @details Markers are assumed independent. For each marker, the possible genotypes
#' for `id` are evaluated and ranked according to how probable they are. It may be wise
#' to try first with `maxPerMarker = 1` to limit computation time, particularly if mutations are modeled.
#'
#' @export
#' @return A list with three components. The first, a data frame with profiles
#' ranked according to the likelihood. The posterior equals the likely for a prior which
#' gives equal probability to all possible profiles. Then follows the most likely profile. Finally,
#' the second to most likely genotypes are given.
#' @author Magnus Dehli Vigeland
#' @examples
#' library(pedprobr)
#' x = nuclearPed(2, father = "FA")
#' m1 = marker(x, "3" = 1, "4" = 1:2, alleles = 1:2, name = "m1")
#' m2 = marker(x, "3" = 1, "4" = 2, alleles = 1:2, name = "m2")
#' m3 = marker(x, "3" = 2, "4" = 2, alleles = 1:2, name = "m3")
#' x = setMarkers(x, list(m1, m2, m3))
#' id = "FA"
#' markers = NULL
#' maxPerMarker = 2
#' rankProfiles(x, id, markers, maxPerMarker)
#' rankProfiles(x, id, markers, maxPerMarker = 1)
#' library(pedmut)
#' mutmod(x, 1:3) = list(model = "equal", rate = 0.002)
#' rankProfiles(x, id, markers, maxPerMarker = Inf)
rankProfiles = function(x, id, markers = NULL, maxPerMarker = Inf, verbose = FALSE) {

  if(!hasMarkers(x))
    stop("No markers attached to pedigree")

  if(is.null(markers))
    markers = seq_len(nMarkers(x))

  # Extract wanted markers
  x = selectMarkers(x, markers)
  nmark = nMarkers(x)

  # Marker names
  mnames = name(x, 1:nmark)
  mnames[is.na(mnames)] = as.character((1:nmark)[is.na(mnames)])

  # Likelihood for each marker
  omdList = lapply(1:nmark, function(i) {
    omd = oneMarkerDistribution(x, ids = id, partialmarker = i, verbose = FALSE)
    a = omd[omd > 0, drop = F]

    if(maxPerMarker < length(a))
      a = sort(a, decreasing = T)[1:maxPerMarker]
    else
      a = sort(a, decreasing = T)[1:length(a)]

    if(!is.array(a)) a = as.array(a)

    # Trick to ensure proper naming in `expand.grid` below
    names(dimnames(a)) = mnames[i]

    if(verbose) {
      print(a)
      cat("\n")
    }

    a
  })
    # Find most likely profile
    maxGenos = lapply(omdList, function(x) x[1])
    maxGenos = unlist(maxGenos)
    names(maxGenos)=paste(mnames, names(maxGenos), sep=":")

    # Find second most likely marginally
    no = unlist(lapply(omdList, length))
    seconds = no > 1
    if (any(seconds)){
      marginal2 = lapply(omdList[seconds], function(x) x[2])
      marginal2 = unlist(marginal2)
      names(marginal2)=paste(mnames[seconds], names(marginal2), sep=":")
    } else
        marginal2 = NA


    # Array with total posterior probs
      pp = Reduce(`%o%`, omdList)

      # Transform dimnames into data frame with possible profiles
      profiles = expand.grid(dimnames(pp))

      # Add probabilities
      prob = as.numeric(pp)
      profiles$likelihood = prob
      # Sort
      profiles = profiles[order(prob, decreasing = T), , drop = F]
  # }
  row.names(profiles) = NULL

  list("rankedProfiles" = profiles, "maxProfile" = maxGenos,
       "no2" = marginal2)
}
