########################################################################################################################
## options.R
## created: 2014-12-04
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Routines for manipulating the echipp options.
########################################################################################################################

## G L O B A L S #######################################################################################################

## Names of all recognized echipp options
ECHIPP.OPTIONS <- c("echipp.aligner", "echipp.bowtie", "echipp.samtools", "echipp.macs",
	"echipp.cluster", "echipp.verbose")

## F U N C T I O N S ###################################################################################################

echipp.load.options <- function(fname = system.file("data/options.RData", package = "echipp")) {
	if (file.exists(fname)) {
		load(fname)
	} else {
		echipp.options <- list()
	}
	echipp.options
}

########################################################################################################################

#' echipp.reload.options
#'
#' Resets all \pkg{echipp} options to their default values, if such are available.
#'
#' @return Invisibly, \code{list} containing the previous values of the reloaded options; an empty list if there are
#'         no default option values stored in the current installation of \pkg{echipp}.
#' @author Yassen Assenov
#' @export
echipp.reload.options <- function() {
	echipp.options <- echipp.load.options()
	if (length(echipp.options) != 0) {
		return(invisible(do.call(options, echipp.options)))
	}
	return(invisible(list()))
}
echipp.reload.options()

########################################################################################################################

#' echipp.save.options
#'
#' Saves the currently set \pkg{echipp} options as default.
#'
#' @return Invisibly, \code{TRUE} if all options were saved as default, \code{FALSE} otherwise.
#' @author Yassen Assenov
#' @export
echipp.save.options <- function() {
	echipp.options <- do.call(options, as.list(ECHIPP.OPTIONS))
	fname <- file.path(system.file(package = "echipp"), "data/options.RData")
#	if (all(sapply(echipp.options, is.null))) {
#		return(invisible((!file.exists(fname)) || suppressWarnings(file.remove(fname))))
#	}
	save(echipp.options, file = fname, compression_level = 9L)
	return(invisible(TRUE))
}

########################################################################################################################

#' echipp.is.verbose
#'
#' Checks if verbose output for \pkg{echipp} functions is enabled.
#'
#' @return \code{TRUE} if the flag \code{"echipp.verbose"} is set to \code{TRUE}, or if it is not set but the option
#'         \code{"verbose"} is set to \code{TRUE}.
#' @author Yassen Assenov
#' @noRd
echipp.is.verbose <- function() {
	result <- getOption("echipp.verbose")
	if (is.null(result)) {
		result <- isTRUE(getOption("verbose"))
	}
	result
}

########################################################################################################################

#' echipp.get.location
#'
#' Reads the location (directory) of a tools specified as an \pkg{echipp} option, or attempts to extract it from the
#' system path.
#'
#' @param option.name        Option name as a one-element \code{character} vector.
#' @param option.description Short description of the option as a one-element \code{character} vector. This is used
#'                           in constructing an error message in case the tool cannot be located.
#' @param executable.name    Optionally, name of the executable file. If this parameter is set, this function validates
#'                           that the specified file exists and is executable.
#' @param attempt.create     Flag indicating if the directory must be created in case it doesn't exist. This parameter
#'                           is ignored when \code{executable.name} is set.
#' @return Location of the tool as a one-element \code{character} vector.
#' @author Yassen Assenov
#' @noRd
echipp.get.location <- function(option.name, option.description, executable.name = NULL, attempt.create = FALSE) {
	result <- getOption(option.name)
	if (is.null(result)) {
		if (!is.null(executable.name)) {
			result <- suppressMessages(suppressWarnings(
				system(paste('which', executable.name), intern = TRUE, show.output.on.console = FALSE)))
			if (length(result) == 0 || isTRUE(attr(result, "status") == 1)) {
				result <- NULL
			}
		}
		if (is.null(result)) {
			txt <- ifelse(is.null(executable.name), "", " and could not be inferred")
			stop(paste0("Missing ", option.description, txt, "; set option ", option.name))
		}
	}
	if (!(is.character(result) && length(result) == 1 && isTRUE(result != ""))) {
		stop(paste0("Invalid value for ", option.description, "; change option ", option.name))
	}
	if (is.null(executable.name)) {
		if ((!file.exists(result)) && (!dir.create(result, FALSE, TRUE))) {
			stop(paste0("Could not create ", option.description, "; change option ", option.name))
		}
	} else if (!isTRUE((file.info(file.path(result, executable.name))[1, "mode"] & "111") != "0")) {
		stop(paste0("Invalid ", option.description, "; update option ", option.name))
	}
	result
}
