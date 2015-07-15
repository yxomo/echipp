#' EchippSet Class
#'
#' Basic class for storing information on ChIP-seq samples.
#'
#' @details
#' This class contains and works on a sample annotation table. The table contains, among other information, links to the
#' files with ChIP-seq reads, peaks or other inferred regions.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{info}}{Sample annotations (phenotypic and processing data) in the form of a \code{data.frame}.}
#'   \item{\code{colors}}{Mappings from categories defined in the sample annotation table to colors.}
#' 	 \item{\code{assembly}}{\code{character} vector of length one, specifying the genome assembly which the dataset is
#' 	 	linked to, e.g. \code{"mm10"}.}
#' }
#'
#' @section Methods and Functions:
#' \describe{
#'   \item{\code{\link{echipp.initialize}}}{Initializes a new dataset from CSV file(s).}
#'   \item{\code{\link[=samples,EchippSet-method]{samples}}}{Gets the identifiers of all samples in the dataset.}
#'   \item{\code{\link[=assembly,EchippSet-method]{assembly}}}{Get the genome assembly associated with the dataset.}
#'   \item{\code{\link{echipp.dataset.files}}}{...}
#'   \item{\code{\link{echip.set.colors}}}{Adds or replaces color mappings to a dataset.}
#' }
#'
#' @name EchippSet-class
#' @rdname EchippSet-class
#' @author Yassen Assenov
#' @exportClass EchippSet
setClass("EchippSet",
	representation(
		info = "data.frame",
		colors = "list",
		genome = "character"),
	prototype(
		info = data.frame(),
		colors = list(),
		genome = as.character(NA)),
	validity = echipp.valid.EchippSet,
	package = "echipp")
