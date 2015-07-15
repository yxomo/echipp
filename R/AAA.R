########################################################################################################################
## AAA.R
## created: 2015-05-01
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Global variables and functions validating commonly used data structures.
########################################################################################################################

## G L O B A L S #######################################################################################################

#' Easy ChIP-seq analysis Pipeline
#'
#' Echipp facilitates the analysis and comparison of ChIP-seq datasets.
#'
#' @docType package
#' @name echipp
NULL

RECOGNIZED.COLUMNS <- c("ID" = "character", "Name" = "character",
	"Original fastq file 1" = "character", "Original fastq file 2" = "character",
	"Processed fastq file 1" = "character", "Processed fastq file 2" = "character", "Fragment length" = "integer",
	"Input" = "character", "SAM file" = "character", "BAM file" = "character", "Track file" = "character",
	"Positive peak file" = "character", "Negative peak file" = "character", "MACS model file" = "character")

FILE.LOCATIONS <- c("Original fastq file 1" = "fastq-original", "Processed fastq file 1" = "fastq-processed",
	"SAM file" = "sam", "BAM file" = "bam", "Track file" = "tracks", "Positive peak file" = "peaks",
	"Negative peak file" = "peaks", "MACS model file" = "macs")
FILE.ENDINGS <- c("fastq.gz", "sam", "bam", "bedGraph.gz", "pos.txt.gz", "neg.txt.gz", "txt")
names(FILE.ENDINGS) <- names(FILE.LOCATIONS)[-1]

########################################################################################################################

## Regular expression for recognizing fastq files
REGEX.FASTQ <- c("^(.*/)?", "/([^/]+)\\.fastq(\\.gz)?$")

## Regular expression for recognizing peak files
REGEX.PEAKS <- c("^(.*/)?", "/([^/]+)\\.(pos|neg)\\.(bed|txt)(\\.gz)?$")

## F U N C T I O N S ###################################################################################################

#' echipp.valid.annotation
#'
#' Validates that the given data frame matches the requirements for a sample annotation table.
#'
#' @param tbl    Table to be validated, given as a non-empty \code{data.frame}.
#' @param modify Flag specifying if column types and missing values could be converted to match the style expected in
#'               \linkS4class{EchippSet} objects.
#' @return Invisibly, the (possibly modified) sample annotation table.
#' @author Yassen Assenov
#' @noRd
echipp.valid.annotation <- function(tbl, modify = FALSE) {
	## Validate all required columns are present
	i <- setdiff(c("ID"), colnames(tbl))
	if (length(i) != 0) {
		return(paste0("Missing required column(s):", paste(i, collapse = ", ")))
	}
	for (cname in intersect(names(RECOGNIZED.COLUMNS), colnames(tbl))) {
		x <- tbl[, cname]
		if (RECOGNIZED.COLUMNS[cname] == "character") {
			x <- as.character(x)
			if (cname == "ID") {
				if (!(isTRUE(all(x != "")) && anyDuplicated(x) == 0)) {
					return("Missing or duplicated IDs")
				}
				if (!all(grepl("^[0-9A-Za-z_()+-]+$", x))) {
					return("Invalid values in column ID")
				}
			} else if (cname == "Name") {
				if (!(isTRUE(all(x != "")) && anyDuplicated(x) == 0)) {
					return("Missing or duplicated sample names")
				}
			} else {
				x[x == ""] <- NA
			}
		} else if (RECOGNIZED.COLUMNS[cname] == "integer") {
			if (all(is.na(x)) || (is.double(x) && isTRUE(all(na.omit(x) == as.integer(na.omit(x)))))) {
				x <- as.integer(x)
			} else if (!is.integer(x)) {
				return(paste("Invalid values in column", cname))
			}
			if (cname == "Fragment length") {
				if (isTRUE(any(x <= 0))) {
					return(paste("Non-positive values in column", cname))
				}
			}
		}
		if (modify) {
			tbl[, cname] <- x
		}
	}

	if (modify) {
		rownames(tbl) <- tbl[, "ID"]
		## Add column Name
		if (!("Name" %in% colnames(tbl))) {
			tbl[, "Name"] <- tbl[, "ID"]
		}
	}

	## Check information on inputs
	if ("Input" %in% colnames(tbl)) {
		i <- which(!is.na(tbl[, "Input"]))
		if (length(i) != 0) {
			i <- i[!(tbl[i, "Input"] %in% rownames(tbl))]
			if (length(i) != 0) {
				i <- paste(tbl[i, "Input"], collapse = ", ")
				return(paste("The following ID(s) in column Input not found:", i))
			}
		}
	}

	## Check information on single-end and paired-end reads
	if (all(c("Original fastq file 1", "Original fastq file 2", "Fragment length") %in% colnames(tbl))) {
		i <- which(is.na(tbl[, "Original fastq file 2"]) & is.na(tbl[, "Fragment length"]))
		if (length(i) != 0) {
			i <- paste(rownames(tbl), collapse = ", ")
			return(paste("The following ID(s) are missing paired ends or fragment lengths:", i))
		}
		i <- which((!is.na(tbl[, "Original fastq file 2"])) & (!is.na(tbl[, "Fragment length"])))
		if (length(i) != 0) {
			txt <- paste(rownames(tbl), collapse = ", ")
			warning(paste("Fragment lengths for the following ID(s) are ignored:", txt))
		}
	}

	invisible(tbl)
}

########################################################################################################################

#' echipp.valid.EchippSet
#'
#' Validation method for EchippSet objects.
#'
#' @param object Object of type \linkS4class{EchippSet} to be tested.
#' @return \code{TRUE} if \code{object} is a valid dataset; one-element \code{character} vector with an error message
#'         otherwise.
#' @author Yassen Assenov
#' @noRd
echipp.valid.EchippSet <- function(object) {
	tbl <- object@info
	if (is.null(tbl)) {
		return("Missing annotation table")
	}
	result <- echipp.valid.annotation(tbl)
	if (is.character(result)) {
		return(result)
	}
	return(TRUE)
}

########################################################################################################################

#' echipp.update.annotation.columns
#'
#' Updates those columns in the given annotation table that relate to the locations of recognized sample-specific files.
#'
#' @param tbl             Sample annotation table passing the validity criteria of
#'                        \code{\link{echipp.valid.annotation}}.
#' @param working.dir     Working directory, i.e. location of the file from which the table was loaded.
#' @param genome.assembly Genome assembly of the dataset.
#' @return The (possibly modified) sample annotation table.
#' @author Yassen Assenov
#' @noRd
echipp.update.annotation.columns <- function(tbl, genome.assembly, working.dir = getwd()) {
	cname.reference <- names(FILE.LOCATIONS[c(1, 2)])[names(FILE.LOCATIONS[c(1, 2)]) %in% colnames(tbl)][1]
	if (is.na(cname.reference)) {
		return(tbl)
	}
	regex.reference <- paste0(REGEX.FASTQ[1], FILE.LOCATIONS[cname.reference], REGEX.FASTQ[2])
	if (!all(grepl(regex.reference, tbl[, cname.reference]))) {
		return(tbl)
	}
	convert.to.absolute <- function(pth) {
		pth <- as.character(pth)
		i <- which(pth == "")
		if (length(i) != 0) {
			pth[i] <- as.character(NA)
		}
		i <- which(!(is.na(pth) | grepl("^(/|[A-Za-z]:)", pth)))
		if (length(i) != 0) {
			pth[i] <- normalizePath(paste0(working.dir, '/', pth[i]), '/', FALSE)
		}
		pth
	}
	for (cn in intersect(c(cname.reference, sub("1$", "2", cname.reference)), colnames(tbl))) {
		tbl[, cn] <- convert.to.absolute(tbl[, cn])
	}
	apply.regex <- function(cn, regex.target, cname = cname.reference) {
		if (echipp.is.verbose()) {
			message(paste0("Added column '", cn, "' based on '", cname, "'"))
		}
		gsub("//", "/", gsub(regex.reference, regex.target, tbl[, cname]), fixed = TRUE)
	}
	for (cn in setdiff(names(FILE.LOCATIONS[-1]), cname.reference)) {
		if (cn %in% colnames(tbl)) {
			tbl[, cn] <- convert.to.absolute(tbl[, cn])
		} else {
			if (cn == "Processed fastq file 1") {
				regex.target <- paste0("\\1/", FILE.LOCATIONS[cn], "/", genome.assembly, "/\\2.", FILE.ENDINGS[cn])
				tbl[[cn]] <- apply.regex(cn, regex.target)
				if ((!("Processed fastq file 2" %in% colnames(tbl))) && ("Original fastq file 2" %in% colnames(tbl))) {
					tbl[["Processed fastq file 2"]] <- apply.regex(cn, regex.target, "Original fastq file 2")
				}
			} else {
				regex.target <- paste0("\\1/", FILE.LOCATIONS[cn], "/", genome.assembly, "/")
				tbl[[cn]] <- paste0(apply.regex(cn, regex.target), rownames(tbl), ".", FILE.ENDINGS[cn])
			}
		}
	}
	return(tbl)
}
