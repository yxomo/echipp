########################################################################################################################
## macs.R
## created: 2014-12-04
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions running MACS and summarizing its output.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' echipp.run.peak.calling
#'
#' ...
#'
#' @param dataset         Dataset of interest as an object of type \linkS4class{EchippSet}.
#' @param job.name.prefix Prefix of the job name given as one-element \code{character} vector. This parameter is used
#'                        only if cluster environment is enabled.
#' @return ...
#'
#' @author Yassen Assenov
#' @export
echipp.run.peak.calling <- function(dataset, job.name.prefix = 'Echipp_peak_calling_') {
	## Validate parameters
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	if(!(is.character(job.name.prefix) && length(job.name.prefix) == 1 && (!is.na(job.name.prefix)))) {
		stop("Invalid value for job.name.prefix")
	}
	dir.logs <- echipp.get.location("echipp.logs", "log directory", NULL, TRUE)
	macs.location <- echipp.get.location("echipp.macs", "MACS location", "macs14")
	genome.assembly <- sub("hg", "hs", substr(assembly(dataset), 1, 2), fixed = TRUE)

	result <- integer()
	for (sample.id in samples(dataset)) {
		control.id <- dataset@info[sample.id, "Input"]
		if (is.na(control.id)) {
			next
		}
		fnames.output <- c(
			echipp.get.files(dataset, sample.id, "peaks"),
			echipp.get.files(dataset, sample.id, "track"),
			echipp.get.files(dataset, sample.id, "macs"))
		if (is.na(fnames.output[1])) {
			warning(paste("Positive peak file not defined for", sample.id))
		} else if (file.exists(fnames.output[1])) {
			if (echipp.is.verbose()) {
				message(paste("Positive peak file for", sample.id, "already present"))
			}
		} else {
			fname.target <- echipp.get.files(dataset, sample.id, "bam")
			fname.control <- echipp.get.files(dataset, control.id, "bam")
			if (is.na(fname.target)) {
				warning(paste("BAM file not defined for", sample.id))
				next
			}
			if (!file.exists(fname.target)) {
				warning(paste("BAM file is missing for", sample.id))
				next
			}
			if (is.na(fname.control)) {
				warning(paste("BAM file not defined for the input of", sample.id))
				next
			}
			if (!file.exists(fname.control)) {
				warning(paste("BAM file is missing for the input of", sample.id))
				next
			}
			if (!all(echipp.create.dir(fnames.output))) {
				warning(paste("Could not create one or more output directories for", sample.id))
				next
			}
			job.name <- paste0(job.name.prefix, sample.id)
			job.size <- as.integer(ceiling(max(file.info(c(fname.target, fname.control))[, "size"]) / 104856700))
			job.options <- c("-l" =
				sprintf("mem=%d00m,walltime=%d:%02d:00", (job.size * 2L), (job.size %/% 4), (job.size %% 4L) * 15L))
			flength <- dataset@info[sample.id, "Fragment length"]
			if (is.null(flength)) {
				stop("Paired-end reads are not yet supported")
			}
			fnames.output[is.na(fnames.output)] <- ""
			sparams <- c('templocation' = echipp.create.temporary.dir(dir.logs), 'macslocation' = macs.location,
				'treatmentfile' = fname.target, 'controlfile' = fname.control,
				'id' = sample.id, 'assembly' = genome.assembly, 'fragmentlength' = flength,
				'pospeakfile' = fnames.output[1], 'negpeakfile' = fnames.output[2],
				'trackfile' = fnames.output[3], 'macsmodelfile' = fnames.output[4])
			result <- c(result, echipp.start.script("macs.sh", sparams, job.name, job.options))
		}
	}

	invisible(echipp.start.scripts(result))
}

########################################################################################################################

echipp.plot.macs.model <- function(fname, y.step = 0.05) {
#fname <- "D:/Datasets/ChIP-seq/E-MTAB-1506/macs/mm9/H3K9ac_sham.txt"
#y.step <- 0.05
	tbl <- tryCatch(read.csv(fname, header = FALSE, quote = ""), error = function(err) { NULL })
	if (is.null(tbl)) {
		stop(paste("Could not load table from", fname))
	}
	if (!identical(unname(sapply(tbl, class)), c("factor", "numeric"))) {
		stop("Unexpected table structure")
	}
	N <- tapply(tbl[, 2], tbl[, 1], length)
	if (!identical(names(N), c("m", "p", "s"))) {
		stop("Unexpected table structure; unknown factors")
	}
	if (!all(N == N[1])) {
		stop("Unexpected table structure; differing number of observations")
	}
	tbl[, 1] <- factor(as.character(tbl[, 1]), levels = c("p", "m", "s"))
	levels(tbl[, 1]) <- c("forward", "reverse", "shifted")
	colors.z <- c("#800000", "#000080", "#000000")
	names(colors.z) <- levels(tbl[, 1])
	N <- (N[1] - 1) / 2L
	tbl[, 3] <- rep(seq.int(-N,N), 3)
	colnames(tbl) <- c("z", "y", "x")
	maxs <- tapply(1:nrow(tbl), tbl[, 1], function(i) { median(tbl[tbl[i, 2] == max(tbl[i, 2]), 3]) })
	ymax <- ceiling(max(tbl[, 2]) / y.step) * y.step

	pp <- ggplot(tbl) + aes_string(x = 'x', y = 'y', color = 'z') + geom_line(lwd = 1.5) +
		labs(x = "Distance to middle (bp)", y = "Percentage", color = "Tags") +
		scale_y_continuous(limits = c(0, ymax), labels = seq(0, ymax, by = y.step), expand = c(0, 0)) +
		scale_colour_manual(values = colors.z)
	for (lv in names(maxs)) {
		pp <- pp + geom_vline(xintercept = maxs[lv], linetype = 2, color = colors.z[lv])
	}
	pp <- pp + theme(plot.margin = grid::unit(0.1 + c(0, 1, 0, 0), "in")) +
		theme(legend.position = c(1, 0.5), legend.justification = c(0, 0.5))
	return(pp)
#pdf("C:/Users/assenov/DKFZ/Meetings/2015-07-10-NGS/poster/images/macs_H3K9ac_sham.pdf", width = 7.2, height = 7.2)
#print(pp)
#dev.off()
}

########################################################################################################################

#' echipp.report.peaks
#' 
#' Creates a report summarizing the outcome of peak calling algorithms.
#' 
#' @param dataset    Dataset of interest as an object of type \linkS4class{EchippSet}.
#' @param dir.report Directory to contain the generated report. The report's HTML file is named \code{"peaks.html"}. If
#'                   the given directory is non-existent, this function attempts to create it. Set this to \code{"."} to
#'                   use the current working directory.
#' @return The generated report, invisibly.
#' 
#' @author Yassen Assenov
#' @export
echipp.report.peaks <- function(dataset, dir.report) {
#dir.report <- "D:/Datasets/ChIP-seq/ENCODE/echipp-reports"
#suppressPackageStartupMessages(library(RnBeads))
	if (!inherits(dataset, "EchippSet")) {
		stop("invalid value for dataset")
	}
	if (!(is.character(dir.report) && length(dir.report) == 1 && isTRUE(dir.report != ""))) {
		stop("invalid value for dir.report")
	}
	echipp.require("RnBeads")

	if (file.exists(dir.report)) {
		if (!file.info(dir.report)[1, "isdir"]) {
			stop("invalid value for dir.report; existing file")
		}
	} else if (!rnb.initialize.reports(dir.report)) {
		stop("could not initialize report directory")
	}
	ggplot2::theme_set(ggplot2::theme_bw())

	## Initialize the report
	report <- tryCatch(createReport(file.path(dir.report, "peaks.html"), "Peak Calling", "Peak Calling", "echipp",
			init.configuration = !file.exists(file.path(dir.report, "configuration"))),
		error = function(err) { NULL })
	if (is.null(report)) {
		stop("could not initialize report in dir.reports")
	}
	sample.names <- echipp.sample.names(dataset)
	txt <- c('This report summarizes the results of running <a ',
		'href="http://liulab.dfci.harvard.edu/MACS/">MACS</a> on some or all of the available ',
		'<b>', length(sample.names), '</b> sample', ifelse(length(sample.names) == 1, '', 's'), ' from the study of ',
		'interest. Information on how to interpret the plots shown here is available in the dedicated <a ',
		'href="http://www.genomebiology.com/content/9/9/R137">publication</a> about the tool and its underlying ',
		'algorithms.')
	report <- rnb.add.section(report, "Introduction", txt)

	files.model <- sapply(sample.names, function(sample.id) { echipp.get.files(dataset, sample.id, "macs") })
	files.exs <- file.exists(files.model)
	summary.models <- table(files.exs, useNA = "ifany")
	summary.models <- data.frame(
		"File status" = c("FALSE" = "missing", "TRUE" = "present", "NA" = "not specified")[names(summary.models)],
		"Cases" = as.integer(summary.models), check.names = FALSE, stringsAsFactors = FALSE)
	txt <- c("This section summarizes the MACS models that were constructed during peak finding. The table below gives",
		"an overview of the available files, if any.")
	report <- rnb.add.section(report, "MACS Models", txt)
	rnb.add.table(report, summary.models, row.names = FALSE)
	if (isTRUE(any(files.exs))) {
		rplots <- list()
		for (i in 1:length(files.model)) {
			if (is.na(files.exs[i])) {
				pp <- rnb.message.plot("Model file is not specified")
			} else if (files.exs[i]) {
				pp <- tryCatch(echipp.plot.macs.model(files.model[i]), error = function(er) {
						txt <- ifelse(grepl("Unexpected table", er$message),
							 "Invalid or unsupported model file", "Error while opening model file")
						rnb.message.plot(txt)
					}
				)
			} else {
				pp <- rnb.message.plot("Model file is not generated")
			}
			rplot <- createReportPlot(sprintf("models_macs_%04i", i), report, width = 7.2, height = 7.2)
			print(pp)
			rplots <- c(rplots, off(rplot))
		}
		txt <- "The following figure visualizes the enrichments of forward and reverse reads."
		rnb.add.paragraph(report, txt)
		txt <- "MACS models of peaks formed by forward and reverse reads."
		setting.names <- list("Sample" = sample.names)
		names(setting.names[[1]]) <- sprintf("%04d", 1:length(setting.names[[1]]))
		report <- rnb.add.figure(report, txt, rplots, setting.names)
		rm(rplots, i, pp, rplot, txt, setting.names)
	}

	off(report)	
}