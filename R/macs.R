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
				'id' = sample.id, 'assembly' = substr(assembly(dataset), 1, 2), 'fragmentlength' = flength,
				'pospeakfile' = fnames.output[1], 'negpeakfile' = fnames.output[2],
				'trackfile' = fnames.output[3], 'macsmodelfile' = fnames.output[4])
			result <- c(result, echipp.start.script("macs.sh", sparams, job.name, job.options))
		}
	}

	invisible(echipp.start.scripts(result))
}

########################################################################################################################

echipp.plot.macs.model <- function(fname, y.step = 0.05) {
	tbl <- tryCatch(read.csv(fname, header = FALSE, quote = ""), error = function(err) { NULL })
	if (is.null(tbl)) {
		stop(paste("Could not load table from", fname))
	}
#	tbl <- read.csv("D:/Datasets/ChIP-seq/E-MTAB-1506/macs/mm9/H3K9ac_sham.txt", header = FALSE, quote = "")
	if (!identical(unname(sapply(tbl, class)), c("factor", "numeric"))) {
		stop("Unexpected table structure")
	}
	N <- tapply(tbl[, 2], tbl[, 1], length)
	if (!identical(names(N), c("m", "p", "s"))) {
		## TODO: Unknown factors
	}
	if (!all(N == N[1])) {
		## TODO: Unequal number of observations
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

	pp <- ggplot(tbl) + aes(x = x, y = y, color = z) + geom_line(lwd = 1.5) +
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
