## F U N C T I O N S ###################################################################################################

#' echipp.require
#'
#' Attempts to load an optional package and fails with an error if loading is not successfully completed.
#'
#' @param package Name of the package to load as a one-element \code{character} vector.
#'
#' @author Yassen Assenov
#' @noRd
echipp.require <- function(package) {
	if (!suppressWarnings(suppressMessages(require(package, quietly = TRUE, character.only = TRUE)))) {
		stop(paste("Missing required library", package))
	}
}

########################################################################################################################

#' echipp.load.assembly
#'
#' Loads the available annotation for the given genome assembly.
#'
#' @param assembly.code Genome assembly of interest
#' @return The genome annotation as a \code{list} of three elements.
#' @author Yassen Assenov
#' @noRd
echipp.load.assembly <- function(assembly.code) {
	fname <- system.file(paste0("data/", assembly.code, ".RDS"), package = "echipp")
	if (!file.exists(fname)) {
		stop("unsupported assembly")
	}
	genome.annotation <- tryCatch(readRDS(fname), error = function(er) { NULL })
	if (is.null(genome.annotation)) {
		stop(paste("INTERNAL ERROR: Could not load", fname))
	}
	genome.annotation
}

########################################################################################################################

#' echipp.load.table
#'
#' Attempts to load a table from the specified file.
#'
#' @param fname ...
#' @param ... ...
#' @return ...
#'
#' @author Yassen Assenov
#' @export
echipp.load.table <- function(fname, ...) {
	if (!(is.character(fname) && length(fname) == 1 && isTRUE(fname != ""))) {
		stop("Invalid value for fname")
	}
	fun.params <- c(list(file = fname), list(...))
	if (grepl("\\.xls(x?)$", tolower(fname))) {
		echipp.require("openxlsx")
		fun.name <- "read.xlsx"
		names(fun.params)[1] <- "xlsxFile"
	} else if (grepl("\\.csv$", tolower(fname))) {
		fun.name <- "read.csv"
	} else if (grepl("\\.txt$", tolower(fname))) {
		fun.name <- "read.txt"
	} else {
		stop("Invalid value for fname; unknown file type")
	}
	tbl <- tryCatch(suppressWarnings(do.call(fun.name, fun.params)), error = function(err) { NULL })
	if (is.null(tbl)) {
		stop(paste("Cannot open file", fname))
	}
	tbl
}

########################################################################################################################

#' echipp.create.dir
#'
#' Creates, if necesessary, the directories in which the given files should be saved.
#'
#' @param fnames Full paths specifying one or more files.
#' @return Invisibly, \code{logical} vector of the same length as \code{fnames}. Every element in this vector is
#'         \code{TRUE} if the directory already exists or was successfully created, and \code{FALSE} otherwise.
#' @author Yassen Assenov
#' @noRd
echipp.create.dir <- function(fnames) {
	invisible(sapply(dirname(fnames), function(dn) {
			ifelse(dn == "." || file.exists(dn), TRUE, dir.create(dn, showWarnings = FALSE, recursive = TRUE))
		}
	))
}

########################################################################################################################

#' echipp.create.temporary.dir
#'
#' Attempts to generate name (and optionally, create) a new temporary directory.
#' 
#' @param path   Full path of a directory to host the new temporary directory.
#' @param create Flag indicating if this function should attempt to create a directory with the newly generated name.
#' @return Full path of the new temporary directory. If \code{create} is \code{TRUE}, this directory is guaranteed to be
#'         be existent and empty; otherwise, this is a non-existent directory.
#' @author Yassen Assenov
#' @noRd
echipp.create.temporary.dir <- function(path, create = TRUE) {
	LTRS <- c(0:9, LETTERS)
	repeat {
		result <- paste0(LTRS[replicate(8L, sample.int(length(LTRS), 1L))], collapse = "")
		result <- file.path(path, result)
		if (!file.exists(result)) {
			break
		}
	}
	if (create) {
		if (!dir.create(result, FALSE, TRUE)) {
			stop("Could not create temporary directory")
		}
	}
	result
}

########################################################################################################################

#' echipp.get.script
#'
#' ...
#'
#' @param fname ...
#' @return ...
#' @author Yassen Assenov
#' @noRd
echipp.get.script <- function(fname) {
	if (Sys.info()["sysname"] == "Linux") {
		fname.full <- paste0("bin/linux_x86.64/", fname)
	} else if (Sys.info()["sysname"] == "Darwin") {
		fname.full <- paste0("bin/macOSX.i386/", fname)
	} else {
		stop("Unsupported operating system")
	}
	tryCatch(system.file(fname.full, package = "echipp", mustWork = TRUE),
		error = function(e) { stop(paste("Internal error in echipp: required script", fname, "not found")) })
}

########################################################################################################################

#' echipp.start.job
#'
#' Starts a \code{qsub} job on a cluster.
#'
#' @param fname       File to be run in the job, usually a shell script file.
#' @param job.name    Job name.
#' @param job.options Optionally, a named \code{character} vector listing additional job options and their values.
#' @param variables   Optionally, a named \code{character} vector containing variable names and values to be passed
#'                    to the script.
#' @return ...
#' @author Yassen Assenov
#' @noRd
echipp.start.job <- function(fname, job.name, job.options = character(), variables = character()) {
	cmd <- c("-N" = job.name, job.options)
	if (length(variables) != 0) {
		cmd <- c(cmd, c("-v" = paste(names(variables), paste0("'", unname(variables), "'"), sep = "=", collapse = ",")))
	}
	cmd <- paste("qsub", paste(names(cmd), unname(cmd), collapse = " "), fname)
	if (echipp.is.verbose()) {
		cat(paste0(cmd, "\n"))
	}
	if (Sys.info()["sysname"] == "windows") {
#		system(cmd, wait = FALSE, show.output.on.console = FALSE)
	} else {
		system(cmd, wait = FALSE)
	}
	Sys.sleep(0.5)
	return(invisible(0L))
}

########################################################################################################################

#' echipp.start.script
#'
#' ...
#'
#' @param script.file    ...
#' @param script.options ...
#' @param job.name       ...
#' @param job.options    ...
#' @param dir.outcome    ...
#' @return ...
#' @author Yassen Assenov
#' @noRd
echipp.start.script <- function(script.file, script.options, job.name, job.options = character(),
		dir.outcome = getOption("echipp.logs")) {
	fname <- echipp.get.script(script.file)
	if (isTRUE(getOption("echipp.cluster") == "qsub")) {
		job.options <- c('-o' = dir.outcome, job.options)
		echipp.start.job(fname, job.name, job.options, script.options)
	} else {
		txt.command <- paste0(names(script.options), '="', script.options, '"', collapse = " ")
		txt.command <- paste0(txt.command, ' "', fname, '"')
	}
}

########################################################################################################################

#' echipp.start.scripts
#'
#' Starts the commands generated by \code{echipp.start.script}.
#'
#' @param \code{vector} of return values from \code{echipp.start.script}.
#' @return An \code{integer} error code (\code{0} for success).
#'
#' @author Yassen Assenov
#' @noRd
echipp.start.scripts <- function(txt) {
	if (is.character(txt)) {
#		result <- system(paste0("#!/bin/bash\n\n", paste(txt, collapse = "\n\n")))
		result <- system(paste(txt, collapse = "\n\n"))
	} else {
		result <- 0L
	}
	result
}

########################################################################################################################

#' echipp.assemblies
#'
#' Gets the supported assemblies.
#'
#' @return List of all supported assemblies in the form of a \code{character} vector.
#'
#' @author Yassen Assenov
#' @export
echipp.assemblies <- function() {
	regex.assembly <- "^(.+)\\.RDS$"
	result <- dir(system.file("data", package = "echipp"), pattern = regex.assembly)
	result <- setdiff(gsub(regex.assembly, "\\1", result), c("band.colors", "options"))
	result
}

########################################################################################################################

#' echipp.chromosomes
#'
#' Gets the supported chromosomes for the given assembly.
#'
#' @param assembly.code Genome assembly of interest.
#' @return List of supported chromosomes as a \code{character} vector.
#'
#' @author Yassen Assenov
#' @export
echipp.chromosomes <- function(assembly.code) {
	if (!(is.character(assembly.code) && length(assembly.code) == 1 && isTRUE(assembly.code != ""))) {
		stop("invalid value for assembly.code")
	}

	names(echipp.load.assembly(assembly.code)$Bands)
}

########################################################################################################################

#' echipp.chromosome.lengths
#'
#' Gets the lengths of the supported chromosomes for the given assembly.
#'
#' @param assembly.code Genome assembly of interest.
#' @return List of chromosome lengths as an \code{integer} vector. The element names of this vector are the chromosomes.
#'
#' @author Yassen Assenov
#' @export
echipp.chromosome.lengths <- function(assembly.code) {
	if (!(is.character(assembly.code) && length(assembly.code) == 1 && isTRUE(assembly.code != ""))) {
		stop("invalid value for assembly.code")
	}

	sapply(echipp.load.assembly(assembly.code)$Bands, function(x) { x[nrow(x), 1] })
}

########################################################################################################################

#' echipp.dataset.files
#'
#' Extracts the dedicated file for a targeted sample in the specified dataset.
#'
#' @param dataset Dataset of interest as an object of type \linkS4class{EchippSet}.
#' @return Invisibly, \code{logical} \code{matrix} in which rows indicate samples and columns - respective files.
#'
#' @author Yassen Assenov
#' @export
echipp.dataset.files <- function(dataset) {
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	rcolumns <- grep(" file", names(RECOGNIZED.COLUMNS), value = TRUE, fixed = TRUE)
	rcolumns <- intersect(rcolumns, colnames(dataset@info))
	result <- file.exists(unlist(dataset@info[, rcolumns], use.names = FALSE))
	result <- matrix(result, nrow = nrow(dataset@info),
		ncol = length(rcolumns), dimnames = list(rownames(dataset@info), rcolumns))
	if (echipp.is.verbose()) {
		mm <- matrix("", nrow = nrow(result) + 1, ncol = ncol(result) + 1L)
		cnames <- c("", sub(" file", "", colnames(result), fixed = TRUE))
		cnames[1] <- paste(rep(" ", max(nchar(rownames(result)))), collapse = "")
		mm[1, ] <- cnames
		mm[-1, 1] <- sprintf(paste0("%", nchar(cnames[1]), "s"), rownames(result))
		for (i in 1:ncol(result)) {
			mm[-1, i + 1] <- sprintf(paste0("%", nchar(cnames[i + 1]), "s"), ifelse(result[, i], "x", "."))
		}
		mm <- apply(mm, 1, paste, collapse = " ")
		invisible(mapply(message, mm, USE.NAMES = FALSE))
	}
	invisible(result)
}

########################################################################################################################

#' echipp.get.files
#'
#' Extracts the dedicated file for a targeted sample in the specified dataset.
#'
#' @param dataset    Dataset of interest as an object of type \linkS4class{EchippSet}.
#' @param sample.id  Targeted sample in the dataset; this must be an element of \code{samples(dataset)}.
#' @param file.type  Request type of file(s); this must be one of \code{"fastq-original"}, \code{"fastq-processed"},
#'                   \code{"fastq"}, \code{"sam"}, \code{"bam"}, \code{"track"}, \code{"peaks"} or \code{"macs"}.
#' @param named      Flag indicating if the returned vector should contain names. If this is set to \code{TRUE}, the
#'                   names of the elements are constructed from \code{sample.id}.
#' @return File name(s), in the form of a \code{character} vector, for the requested file types of the targeted sample.
#'
#' @author Yassen Assenov
#' @export
echipp.get.files <- function(dataset, sample.id, file.type = "fastq", named = FALSE) {
	## Validate parameters
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	if (!(is.character(sample.id) && length(sample.id) == 1 && isTRUE(sample.id != ""))) {
		stop("Invalid value for sample.id")
	}
	if (!(sample.id %in% samples(dataset))) {
		stop(paste("Unknown sample", sample.id))
	}
	if (!(is.character(file.type) && length(file.type) == 1 && isTRUE(file.type != ""))) {
		stop("Invalid value for file.type")
	}
	if (!(is.logical(named) && length(named) == 1 && (!is.na(named)))) {
		stop("Invalid value for named")
	}

	tbl.value <- function(cname) {
		if (cname %in% colnames(dataset@info)) {
			dataset@info[sample.id, cname]
		} else {
			as.character(NA)
		}
	}
	if (file.type == "fastq") {
		result <- c(tbl.value("Original fastq file 1"), tbl.value("Original fastq file 2"),
			tbl.value("Processed fastq file 1"), tbl.value("Processed fastq file 2"))
		if (is.na(result[2])) {
			if (is.na(result[4])) {
				result <- result[c(1, 3)]
			} else {
				result <- result[-2]
			}
		} else if (is.na(result[4])) {
			result <- result[-4]
		}
	} else if (file.type == "fastq-original") {
		result <- c(tbl.value("Original fastq file 1"), tbl.value("Original fastq file 2"))
		if (is.na(result[2])) { result <- result[1] }
	} else if (file.type == "fastq-processed") {
		result <- c(tbl.value("Processed fastq file 1"), tbl.value("Processed fastq file 2"))
		if (is.na(result[2])) { result <- result[1] }
	} else if (file.type == "sam") {
		result <- tbl.value("SAM file")
	} else if (file.type == "bam") {
		result <- tbl.value("BAM file")
	} else if (file.type == "track") {
		result <- tbl.value("Track file")
	} else if (file.type == "peaks") {
		result <- c(tbl.value("Positive peak file"), tbl.value("Negative peak file"))
	} else if (file.type == "macs") {
		result <- tbl.value("MACS model file")
	} else {
		stop("Unsupported file.type")
	}
	if (named) {
		if (length(result) == 1) {
			names(result) <- sample.id
		} else {
			names(result) <- paste0(1:length(result), '_', sample.id)
		}
	}
	result
}

########################################################################################################################

#' echipp.sample.names
#' 
#' Gets the sample names in a dataset.
#' 
#' @param dataset         Dataset of interest as an object of type \linkS4class{EchippSet}.
#' @return Sample names as a \code{character} vector.
#' 
#' @author Yassen Assenov
#' @export
echipp.sample.names <- function(dataset) {
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	if ("Name" %in% colnames(dataset@info)) {
		return(dataset@info[, "Name"])
	}
	return(rownames(dataset@info))
}

########################################################################################################################

#' echipp.run.alignment
#'
#' ...
#'
#' @param dataset         Dataset of interest as an object of type \linkS4class{EchippSet}.
#' @param job.name.prefix Prefix of the job name given as one-element \code{character} vector. This parameter is used
#'                        only if cluster environment is enabled.
#' @param keep.sam        Flag indicating if the intermediate SAM files should be kept. If this is set to \code{FALSE}
#'                        (default), all SAM files that are created by this process are removed after convertion to
#'                        BAM.
#' @return ...
#'
#' @author Yassen Assenov
#' @export
echipp.run.alignment <- function(dataset, job.name.prefix = 'Echipp_alignment_', keep.sam = FALSE) {
	## Validate parameters
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	if(!(is.character(job.name.prefix) && length(job.name.prefix) == 1 && (!is.na(job.name.prefix)))) {
		stop("Invalid value for job.name.prefix")
	}
	if (!(is.logical(keep.sam) && length(keep.sam) == 1 && (!is.na(keep.sam)))) {
		stop("Invalid value for keep.sam")
	}
	invisible(echipp.get.location("echipp.logs", "log directory", NULL, TRUE))
	bowtie.location <- echipp.get.location("echipp.bowtie", "bowtie location", "bowtie")
	samtools.location <- echipp.get.location("echipp.samtools", "samtools location", "samtools")

	result <- integer()
	for (sample.id in samples(dataset)) {
		fname.bam <- echipp.get.files(dataset, sample.id, 'bam')
		if (is.na(fname.bam)) {
			warning(paste("BAM file not defined for", sample.id))
		} else if (file.exists(fname.bam)) {
			if (echipp.is.verbose()) {
				message(paste("BAM file for", sample.id, "already present"))
			}
		} else {
			fnames.fastq <- echipp.get.files(dataset, sample.id, 'fastq-processed')
			fname.sam <- echipp.get.files(dataset, sample.id, 'sam')
			if (all(is.na(fnames.fastq))) {
				warning(paste("Processed fastq file not defined for", sample.id))
				next
			}
			if (is.na(fname.sam)) {
				warning(paste("SAM file not defined for", sample.id))
				next
			}
			if (!all(echipp.create.dir(c(fname.sam, fname.bam)))) {
				warning(paste("Could not create BAM file or SAM file directory for", sample.id))
				next
			}
			job.name <- paste0(job.name.prefix, sample.id)
			if (length(fnames.fastq) == 1) {
				bparams <- c('bowtielocation' = bowtie.location, 'samtools' = samtools.location,
					'genome' = assembly(dataset), 'fastqfile' = fnames.fastq, 'bamfile' = fname.bam,
					'samfile' = fname.sam, 'keepsam' = ifelse(keep.sam, 'yes', 'no'))
				result <- c(result, echipp.start.script("bowtie.sh", bparams, job.name))
			} else {
				stop("Paired-end reads are not yet supported")
			}
		}
	}

	invisible(echipp.start.scripts(result))
}

########################################################################################################################

#' echipp.run.fastqc
#'
#' ...
#'
#' @param dataset         Dataset of interest as an object of type \linkS4class{EchippSet}.
#' @param output          Optionally, base directory in which to store the FastQC reports. If this is set to \code{NULL}
#'                        (default), every report is saved in the directory of the corresponding file. Otherwise, the
#'                        reports are saved in the following subdirectories of \code{output}: \code{"fastqc-original"},
#'                        \code{"fastqc-processed"}, \code{"bam"}. Any existing reports will be overwritten.
#' @param job.name.prefix Prefix, given as a one-element \code{character} vector, to use for the qsub job that is
#'                        started. This parameter is used only if clustering environment is enabled.
#' @return ...
#'
#' @author Yassen Assenov
#' @export
echipp.run.fastqc <- function(dataset, output = NULL, job.name.prefix = 'Echipp_fastqc_') {
	## Validate parameters
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	if (!is.null(output)) {
		if (!(is.character(output) && length(output) == 1 && isTRUE(output != ""))) {
			stop("Invalid value for output")
		}
		## Create the output directory if necessary
		if (!echipp.create.dir(file.path(output, "a"))) {
			stop("Could not create output directory")
		}
	}
	if (!(is.character(job.name.prefix) && length(job.name.prefix) == 1 && (!is.na(job.name.prefix)))) {
		stop("Invalid value for job.name.prefix")
	}
	invisible(echipp.get.location("echipp.logs", "log directory", NULL, TRUE))
	fastqc.location <- echipp.get.location("echipp.fastqc", "fastqc location", "fastqc")

	echipp.fastqc.file <- function(fastqc.location, fnames, ftype, job.name.prefix, rdir) {
		for (i in 1:length(fnames)) {
			fname <- unname(fnames[i])
			if (!is.na(fname)) {
				sample.id <- names(fnames)[i]
				if (!file.exists(rdir[i])) {
					if (!dir.create(rdir[i], showWarnings = FALSE, recursive = TRUE)) {
						warning(paste("Could not create output directory", rdir))
						break
					}
				}
				job.name <- paste0(job.name.prefix, ftype, '_', sample.id)
				fparams <- c('fastqclocation' = fastqc.location, 'inputfile' = fname, 'sample' = sample.id,
					'reportdir' = rdir[i])
				echipp.start.script("fastqc.sh", fparams, job.name)
				if (echipp.is.verbose()) {
					message(paste('Started job for', ftype, 'file of', sample.id))
				}
				Sys.sleep(0.5)
			}
		}
	}

	result <- integer()
	for (sample.id in samples(dataset)) {
		## Run FastQC on the fastq file(s)
		for (ftype in paste0('fastq-', c('original', 'processed'))) {
			fnames <- echipp.get.files(dataset, sample.id, ftype, TRUE)
			if (is.null(output)) {
				rdir <- dirname(fnames)
			} else {
				rdir <- rep(file.path(output, ftype), length(fnames))
			}
			result <- c(result, echipp.fastqc.file(fastqc.location, fnames, ftype, job.name.prefix, rdir))
		}

		## Run FastQC on the BAM file
		fname <- echipp.get.files(dataset, sample.id, 'bam', TRUE)
		if (isTRUE(file.exists(fname))) {
			rdir <- file.path(output, 'bam', dataset@genome)
#			result <- c(result, echipp.fastqc.file(fastqc.location, fname, 'bam', job.name.prefix, rdir))
		}
	}

	invisible(echipp.start.scripts(result))
}
