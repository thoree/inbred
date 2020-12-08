#' Blind search table
#'
#' Calculate table ofLRs and posteriors for PO, SI, SI, FC and UN for all pairs
#'
#'
#' @param x A `ped` object or a list of such.
#' @param ids Either a vector with ID labels, or a data frame/matrix with two
#'   columns, where each row contains the ID labels of two individuals. The entries
#'   are coerced to characters, and must match uniquely against the ID labels of
#'   `x`. If `ids` is a vector, it is converted to a matrix containing all pairs.
#'   By default, all individuals of `x` are included.
#' @param markers A numeric indicating which marker(s) to include. If NULL
#'   (default), all markers are used.
#' @return A list with two elements, one data frame  for LRs, one for posteriors
#' with flat prior. The  data frame with 6 columns: `ID1`, `ID2`, `N` (the number of markers
#' with no missing alleles), `andLRs or posteriors.
#' @author Thore Egeland
#' @details A uniform priorisused so far
#'
#' @references
#'
#' To appear
#'
#' @examples
#'
#' ### Example 1
#' library(forrel)
#' x = nuclearPed(2)
#' # Simulate 200 equifrequent SNPs
#' x = markerSim(x, N = 200, alleles = 1:2, verbose = FALSE)
#'
#' # Parent offspring against unrelated
#' tab1 = allPairsTable(x)
#'
#' @importFrom utils combn
#' @export
allPairsTable = function(x, ids = NULL, markers = NULL) {
    if(is.ped(x))
        x = list(x)
    else if(!is.pedList(x))
        stop2("The first argument must be a `ped` object or a list of such")

    if(is.null(markers))
        markers = seq_len(nMarkers(x))

    if(is.null(ids))
        ids = unlist(labels(x))
    if(length(ids) < 2)
        stop2("`ids` must be either a vector of length at least 2, or a matrix/data.frame with 2 columns")

    if(is.vector(ids))
        ids = t.default(combn(ids, 2))
    else if(is.data.frame(ids))
        ids = as.matrix(ids)

    if(!is.matrix(ids) && ncol(ids) == 2)
        stop2("`ids` must be either a vector of length at least 2, or a matrix/data.frame with 2 columns")

    ids_df_list = lapply(seq_len(nrow(ids)), function(i) pedlistMembership(x, ids[i,]))


    res = lapply(ids_df_list, function(ids_df) {
        A = IBDest_getAlleleData(x, ids_df, markers)
        # Remove markers with missing alleles
        if(any(miss <- apply(A, 2, anyNA)))
            A = A[, !miss, drop = FALSE]

        # Likelihood function
        a=A[1,]; b=A[2,]; cc=A[3,]; d=A[4,]; pa=A[5,]; pb=A[6,]; pc=A[7,]; pd=A[8,]
        loglikUN = sum(log(.IBDlikelihood(c(1,0),a,b,cc,d,pa,pb,pc,pd)))
        likUN = exp(loglikUN)
        loglikPO = sum(log(.IBDlikelihood(c(0,0),a,b,cc,d,pa,pb,pc,pd)))
        loglikSI = sum(log(.IBDlikelihood(c(0.25,0.25),a,b,cc,d,pa,pb,pc,pd)))
        loglikHS = sum(log(.IBDlikelihood(c(0.50,0.0),a,b,cc,d,pa,pb,pc,pd)))
        loglikFC = sum(log(.IBDlikelihood(c(0.75,0.0),a,b,cc,d,pa,pb,pc,pd)))
        data.frame(ID1 = ids_df$id[1],
                   ID2 = ids_df$id[2],
                   N = ncol(A),
                   LR.PO = exp(loglikPO)/likUN,
                   LR.SI = exp(loglikSI)/likUN,
                   LR.HS = exp(loglikHS)/likUN,
                   LR.FC = exp(loglikFC)/likUN,
                   LR.UN = 1,
                   stringsAsFactors = FALSE)
    })

    LRs = do.call(rbind, res)
    posteriors = data.frame(LRs[,1:3],
             t(apply(LRs[,-(1:3)], 1, function(x) x/sum(x))))
    colnames(posteriors)[-(1:3)] = c("Pr.PO", "Pr.SI",  "Pr.HS",
                                     "Pr.FC", "Pr.UN")
    list(LRs = LRs, posteriors = posteriors)
}

# Match ID labels against a list of pedigrees
pedlistMembership = function(pedl, ids) {
    ped_match = vapply(pedl,
                       function(p) ids %in% labels(p),
                       FUN.VALUE = logical(length(ids)))

    # Convert to matrix if length(ids == 1)
    ped_match = rbind(ped_match, deparse.level = 0)

    if(any(nomatch <- rowSums(ped_match) == 0))
        stop2("IDs not found in pedigrees: ", ids[nomatch])
    if(any(multimatch <- rowSums(ped_match) > 1))
        stop2("IDs matching multiple pedigrees: ", ids[multimatch])

    pednr = apply(ped_match, 1, which)

    data.frame(pednr = pednr, id = ids, stringsAsFactors = FALSE)
}


IBDest_getAlleleData = function(x, ped_id_df, markers = NULL) {
    stopifnot(is.pedList(x), is.data.frame(ped_id_df))
    pednr1 = ped_id_df$pednr[1]
    pednr2 = ped_id_df$pednr[2]

    # Match IDs with orig.ids.
    id1_int = internalID(x[[pednr1]], ped_id_df$id[1])
    id2_int = internalID(x[[pednr2]], ped_id_df$id[2])

    # Collect alleles and frequencies in a matrix with 8 rows
    # (first 4 alleles a,b,c,d followed by their freqs), one column per marker
    # Note that alleles are internal integers

    if(pednr1 == pednr2) {
        ped = x[[pednr1]]
        A = vapply(ped$MARKERS[markers], function(m) {
            als = c(m[id1_int,], m[id2_int,])
            frq = rep_len(NA_real_, length(als))
            frq[als > 0] = attr(m, 'afreq')[als] # works, since 0's in als are ignored when indexing
            c(als, frq)
        }, FUN.VALUE = numeric(8))
    }
    else {
        ped1 = x[[pednr1]]
        ped2 = x[[pednr2]]
        A = vapply(markers, function(i) {
            m1 = ped1$MARKERS[[i]]
            m2 = ped2$MARKERS[[i]]
            als = c(m1[id1_int,], m2[id2_int,])
            frq = rep_len(NA_real_, length(als))
            frq[als > 0] = attr(m1, 'afreq')[als] # works, since 0's in als are ignored when indexing
            c(als, frq)
        }, FUN.VALUE = numeric(8))
    }

    return(A)
}


.IBDlikelihood = function(k, a, b, cc, d, pa, pb, pc, pd) {
    ### Vectorized function for computing kappa likelihoods, given genotypes for two related individuals
    # k: numeric of length 2 = (kappa0, kappa2)
    # a: vector of positive integers (allele 1 of individual 1)
    # b: vector of positive integers (allele 2 of individual 1)
    # cc: vector of positive integers (allele 1 of individual 2) (avoid overloading of 'c')
    # d: vector of positive integers (allele 2 of individual 2)
    # pa, pb, pc, pd: numeric vectors with frequencies of the above alleles.
    #
    # Note that all input vectors (except k) have the same length = #markers.
    homoz1 = a == b
    homoz2 = cc == d
    mac = a == cc
    mbc = b == cc
    mad = a == d
    mbd = b == d
    g1.fr = 2^(!homoz1) * pa * pb
    g2.fr = 2^(!homoz2) * pc * pd

    # Prob(g1, g2 | unrelated)
    UN = g1.fr * g2.fr

    # Prob(g1, g2 | parent-offspring)
    PO = (.5)^(homoz1+homoz2)*pa*pb*(pd*(mac+mbc) + pc*(mad+mbd))

    # Prob(g1, g2 | monozygotic twins)
    MZ = g1.fr * ((mac & mbd) | (mad & mbc))

    # return likelihoods (Thompson)
    k[1] * UN + (1-k[1]-k[2]) * PO + k[2] * MZ
}
