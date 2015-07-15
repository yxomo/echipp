########################################################################################################################
## colors.R
## created: 2014-12-01
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Defines functions for loading and validating color mappings.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' echipp.fix.colors
#'
#' Validates and standardizes the given color names.
#' 
#' @param x Color names in the form of a \code{character} vector.
#' @return The possibly modified color names.
#' @details This function validates that every elements of \code{x} is either a recognizable color word, as returned by
#'          \code{\link[grDevices]{colors}}, or a solid color value in the form \code{"#XXXXXX"}, where \code{X} is a
#'          hexadecimal symbol.
#' @author Yassen Assenov
#' @noRd
echipp.fix.colors <- function(x) {
	i <- grepl("^#[0-9A-Fa-f]{6}$", x)
	x[i] <- toupper(x[i])
	x[!i] <- tolower(x[!i])
	i <- which(!i)[which(!(x[!i] %in% colors()))]
	if (length(i) != 0) {
		stop(x[i[1]])
	}
	x
}

########################################################################################################################

#' echipp.load.colors
#'
#' Loads color mappings from a 3-column comma-separated value file and transformes them to a list of named vectors.
#'
#' @param fname File storing color mappings in a three-column table with the following column names: \code{Trait},
#'              \code{Category}, \code{Color}. If empty (missing) values are found, this function throws an error.
#' @return Loaded color mapping as a \code{list} of named \code{character} vectors. Each element of the list corresponds
#'         to one trait. The names in the vectors are categories, and the values - respective colors.
#'
#' @author Yassen Assenov
#' @noRd
echipp.load.colors <- function(fname) {
	## Load the table
	tbl <- tryCatch(
		suppressWarnings(read.csv(fname, na.strings = character(), check.names = FALSE, stringsAsFactors = FALSE)),
		error = function(err) { NULL })
	if (is.null(tbl)) {
		stop(paste("Cannot open file", fname))
	}
	
	## Validate table structure
	if (!(setequal(colnames(tbl), c("Trait", "Category", "Color")))) {
		stop("Invalid color table; unexpected structure")
	}
	set.to.character <- function(x) {
		if (!is.character(x)) {
			x <- as.character(x)
		}
		x
	}
	tbl$Trait <- set.to.character(tbl$Trait)
	tbl$Category <- set.to.character(tbl$Category)
	tbl$Color <- set.to.character(tbl$Color)
	if (any(is.na(tbl$Trait) | is.na(tbl$Category) | is.na(tbl$Color))) {
		stop("Invalid color table; missing values")
	}
	if (any(tbl$Trait == "" | tbl$Category == "" | tbl$Color == "")) {
		stop("Invalid color table; missing values")
	}
	
	## Validate and convert color values
	tbl$Color <- tryCatch(echipp.fix.colors(tbl$Color), error = function(err) {
		stop(paste("Invalid color table; unknown color", err$message))
	}
	)
	
	## Transform the colors to list
	tapply(1:nrow(tbl), tbl$Trait, function(i) {
		result <- tbl[i, "Color"]
		names(result) <- tbl[i, "Category"]
		result
	},
	simplify = FALSE)
}

########################################################################################################################

#' echipp.filter.colors
#' 
#' Filters the provided color mappings to traits available in the sample annotation table.
#' 
#' @param color.mappings Mappings from trait categories to colors. This must be a named \code{list} of named
#'                       \code{character} vectors storing the mappings.
#' @param info.data      Sample annotation table of a dataset.
#' @return Filtered color mappings containing only the traits in \code{info.table} for which all categories are mapped
#'         to colors.
#' @author Yassen Assenov
#' @noRd
echipp.filter.colors <- function(color.mappings, info.table) {
	result <- list()
	for (cname in intersect(colnames(info.table), names(color.mappings))) {
		categories.expected <- as.character(na.omit(unique(info.table[, cname])))
		categories.missing <- setdiff(categories.expected, names(color.mappings[[cname]]))
		if (length(categories.missing) == 0) {
			result[[cname]] <- color.mappings[[cname]]
		} else {
			txt <- paste0("Missing color for the following categor",
				ifelse(length(categories.missing) == 1, "y", "ies"), " in trait ", cname, ": ",
				paste(categories.missing, collapse = ", "))
			warning(txt)
		}
	}
	if (echipp.is.verbose()) {
		for (cname in setdiff(colnames(info.table), c(names(RECOGNIZED.COLUMNS), names(color.mappings)))) {
			message(paste0("Trait '", cname, "' not found in the color mapping table"))
		}
	}
	result
}

########################################################################################################################

#' echipp.set.colors
#'
#' Adds or replaces color mappings to the specified Echipp dataset.
#' 
#' @param dataset        Dataset to be updated.
#' @param color.mappings Mappings from trait categories to colors. This must be a named \code{list} of named
#'                       \code{character} vectors storing the mappings.
#' @param clean          Flag indicating if all previous color mappings of \code{dataset} should be removed.
#' @return The modificed dataset.
#' 
#' @author Yassen Assenov
#' @export
echipp.set.colors <- function(dataset, color.mappings, clean = FALSE) {
	## Validate parameters
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	if (!is.list(color.mappings)) {
		stop("Invalid value for color.mappings; expected list")
	}
	is.well.named <- function(x) { (!is.null(names(x))) && isTRUE(all(names(x) != "")) }
	if (length(color.mappings) != 0) {
		if (!is.well.named(color.mappings)) {
			stop("Invalid value for color.mappings; missing trait names")
		}
		if (!all(sapply(color.mappings, function(x) { is.vector(x) && is.character(x) && length(x) != 0 }))) {
			stop("Invalid value for color.mappings; expected non-empty character vector(s)")
		}
		if (!all(sapply(color.mappings, is.well.named))) {
			stop("Invalid value for color.mappings; expected named vector(s)")
		}
		color.mappings <- tryCatch(lapply(color.mappings, echipp.fix.colors), error = function(err) { err$message })
		if (is.character(color.mappings)) {
			stop(paste("Invalid value for color.mappings; unknown color", color.mappings))
		}
	}
	if (!(is.logical(clean) && length(clean) == 1 && (!is.na(clean)))) {
		stop("Invalid value for clean")
	}
	rm(is.well.named, fix.colors)

	## Add or replace all color mappings
	if (clean) {
		dataset@colors <- list()
	}
	color.mappings <- echipp.filter.colors(color.mappings, dataset@info)
	dataset@colors[names(color.mappings)] <- color.mappings
	dataset
}

########################################################################################################################

#' echipp.color.mappings
#' 
#' Extract the mapped and not mapped traits in the annotation table of the given dataset.
#' 
#' @param dataset         Dataset of interest as an object of type \linkS4class{EchippSet}.
#' @return Traits, along with their mapping to colors as a \code{list} of named \code{character} vectors. Furthermore,
#'         the returned object contains attributes \code{"unmappable"} and \code{"unmapped"}, storing the names of all
#'         properties that are inherently not storing categories (not mappable to colors), and not (fully) considered in
#'         the dataset's color mapping tables.
#' 
#' @author Yassen Assenov
#' @export
echipp.color.mappings <- function(dataset) {
	if (!inherits(dataset, "EchippSet")) {
		stop("invalid value for dataset")
	}
	traits <- setdiff(colnames(dataset@info), names(RECOGNIZED.COLUMNS))
	traits <- lapply(dataset@info[, traits], function(x) {
			if (class(x) == "Date") {
				x <- as.character(x)
			}
			if (is.character(x)) {
				return(as.character(unique(na.omit(x))))
			} else if (is.factor(x)) {
				return(levels(x))
			}
			return(NULL)
		}
	)
	## Define unmappable traits
	i <- which(as.logical(sapply(traits, is.null)))
	if (length(i) != 0) {
		traits.unmappable <- names(traits)[i]
		traits <- traits[-i]
	} else {
		traits.unmappable <- character()
	}
	## Define unmapped traits
	cn <- names(dataset@colors)
	if (is.null(cn)) { cn <- character() }
	traits.unmapped <- setdiff(names(traits), cn)
	traits <- traits[intersect(names(traits), cn)]
	for (tname in names(traits)) {
		traits[[tname]] <- dataset@colors[[tname]][intersect(names(dataset@colors[[tname]]), traits[[tname]])]
	}
	attr(traits, "unmappable") <- traits.unmappable
	attr(traits, "unmapped") <- traits.unmapped
	traits
}
