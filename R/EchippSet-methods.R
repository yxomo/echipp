## A C C E S S O R S ###################################################################################################

if (!isGeneric("assembly")) setGeneric("assembly", function(object) standardGeneric("assembly"))

#' assembly-methods
#'
#' Gets the genome assembly.
#'
#' @param object Dataset of interest.
#' @return Genome assembly to which the dataset is linked.
#'
#' @rdname assembly-methods
#' @docType methods
#' @aliases assembly
#' @aliases assembly,EchippSet-method
#' @export
setMethod("assembly", signature(object = "EchippSet"), function(object) { object@genome })

########################################################################################################################

if (!isGeneric("samples")) setGeneric("samples", function(object) standardGeneric("samples"))

#' samples-methods
#'
#' Gets the sample identifiers.
#'
#' @param object Dataset of interest.
#' @return Identifiers of all samples in the dataset, in the form of a \code{character} vector.
#'
#' @rdname samples-methods
#' @docType methods
#' @aliases samples
#' @aliases samples,EchippSet-method
#' @export
setMethod("samples", signature(object = "EchippSet"), function(object) { rownames(object@info) })

########################################################################################################################

setMethod("show", "EchippSet", function(object) {
		cat("Object of class ", class(object), "\n", sep = "")
		tbl <- object@info
		cat(sprintf("%7d samples; using genome %s\n\n", nrow(tbl), object@genome))
		stats <- NULL
		rcolumns <- grep(" file", names(RECOGNIZED.COLUMNS), value = TRUE, fixed = TRUE)
		for (rcolumn in intersect(rcolumns, colnames(tbl))) {
			fnames <- tbl[, rcolumn]
			fnames[!is.na(fnames)]
			if (length(fnames) != 0) {
				fnames <- file.exists(fnames)
				stats <- rbind(stats,
					data.frame(Li = length(fnames), Pr = sum(fnames), Ty = rcolumn, stringsAsFactors = FALSE))
			}
		}
		if (!is.null(stats)) {
			cat(" Listed Present Type\n--------------------\n")
			for (i in 1:nrow(stats)) {
				cat(sprintf("%7d %7d %s\n", stats[i, 1], stats[i, 2], stats[i, 3]))
			}
			cat("\n")
		}
		traits <- echipp.color.mappings(object)
		t1 <- names(traits)
		t2 <- attr(traits, "unmapped")
		if (length(t1) == 0 && length(t2) == 0) {
			cat("There are no mappable traits in the annotation table.\n")
		} else {
			if (length(t1) != 0) {
				cat("Traits mapped to colors: ", paste(t1, collapse = ", "), "\n", sep = "")
			}
			if (length(t2) != 0) {
				cat("Traits not mapped to colors: ", paste(t2, collapse = ", "), "\n", sep = "")
			}
		}
	}
)

## F U N C T I O N S ###################################################################################################

#' echipp.initialize
#'
#' Initializes an \linkS4class{EchippSet} dataset by loading the sample annotation table and, optionally, the color
#' mappings from CSV file(s).
#'
#' @param fname.data    Name of a comma-separated value (csv) file listing all samples and their annotation.
#' @param assembly.code Genome assembly against which the reads are/should be aligned.
#' @param fname.colors  Name of a comma-separated value (csv) file listing all color mappings for categories defined in
#'                      the sample annotation table.
#' @param infer.files   Flag indicating if the table should be enriched with file name columns.
#' @return Newly initialized dataset as an \linkS4class{EchippSet} object.
#'
#' @author Yassen Assenov
#' @export
echipp.initialize <- function(fname.data, assembly.code, fname.colors = NULL, infer.files = FALSE) {
	## Validate the parameters
	if (!(is.character(fname.data) && length(fname.data) == 1 && isTRUE(fname.data != ""))) {
		stop("Invalid value for fname.data")
	}
	if (!(is.character(assembly.code) && length(assembly.code) == 1 && isTRUE(assembly.code != ""))) {
		stop("Invalid value for assembly.code")
	}
	if (!is.null(fname.colors)) {
		if (!(is.character(fname.colors) && length(fname.colors) == 1 && isTRUE(fname.colors != ""))) {
			stop("Invalid value for fname.colors")
		}
	}
	if (!(is.logical(infer.files) && length(infer.files) == 1 && (!is.na(infer.files)))) {
		stop("Invalid value for infer.files")
	}

	## Load the sample annotation table
	tbl <- echipp.load.table(fname.data, check.names = FALSE, stringsAsFactors = FALSE)
	if (echipp.is.verbose()) {
		message(paste("Loaded", nrow(tbl), "records from", fname.data))
	}

	## Validate and enrich the sample annotation table
	tbl <- echipp.valid.annotation(tbl, TRUE)
	if (is.character(tbl)) {
		stop(tbl)
	}
	if (infer.files) {
		tbl <- echipp.update.annotation.columns(tbl, assembly.code, normalizePath(dirname(fname.data)))
	}

	## Load colors
	if (is.null(fname.colors)) {
		color.mappings <- list()
	} else {
		color.mappings <- echipp.load.colors(fname.colors)
		color.mappings <- echipp.filter.colors(color.mappings, tbl)
	}

	if (!(assembly.code %in% echipp.assemblies())) {
		warning("unsupported assembly")
	}
	new("EchippSet", info = tbl, colors = color.mappings, genome = assembly.code)
}
