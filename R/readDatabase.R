#' Reads database of allele frequencies
#'
#' @param character filename
#' @details An tab separated asci file is read and converted to
#' `ped` format. The first line gives the name of of the marker,
#' the next lines the allele name and frequency. The markers are
#' separated by blank lines.
#' @export
#' @examples
#' \dontrun{
#' library(pedtools)
#' db = readDatabase("markerDB.txt")
#' # Attach to a pedigree
#' x = setMarkers(singleton(1), locusAttributes = db)
#' # Specify X chromosome
#' chrom(x, 1:length(db)) = "X"
#' x
#' }

readDatabase = function(filename) {

  # Read as table with two columns
  raw = read.table(filename, sep = "\t", colClasses = c("character", "numeric"),
                   header = FALSE, fill = TRUE, blank.lines.skip = FALSE)

  # First/last lines numbers for each marker
  newMarker = c(1, which(raw$V1 == "") + 1)
  stops = c(newMarker[-1] -  2, nrow(raw))

  # Convert to list of frequency verctors
  res = lapply(seq_along(newMarker), function(i) {
    m = raw[newMarker[i]:stops[i], ]
    afr = setNames(m[-1, 2], m[-1, 1])
  })

  # Add marker names
  names(res) = raw[newMarker, 1]

  # Return
  res
}

