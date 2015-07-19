########################################################################################################################
## fastqc.R
## created: 2014-12-01
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions related to processing the generated FastQC reports on samples.
########################################################################################################################

## G L O B A L S #######################################################################################################

PLOTS.FASTQC <- c(
	"per_base_quality" = "per base sequence quality",
	"per_sequence_quality" = "per sequence quality",
	"per_base_sequence_content" = "per base sequence content",
	"per_base_gc_content" = "per base GC content",
	"per_sequence_gc_content" = "per sequence GC content",
	"per_base_n_content" = "per base N content",
	"sequence_length_distribution" = "sequence length distribution",
	"duplication_levels" = "duplication levels",
	"kmer_profiles" = "K-mer profiles")

SUMMARY.STATES <- c(
	"Basic Statistics",
	"Per base sequence quality",
	"Per sequence quality scores",
	"Per base sequence content",
	"Per base GC content",
	"Per sequence GC content",
	"Per base N content",
	"Sequence Length Distribution",
	"Sequence Duplication Levels",
	"Overrepresented sequences",
	"Kmer Content")

## F U N C T I O N S ###################################################################################################

echipp.fastqc.summarize <- function(dataset, dir.report, dir.fastq.original = NULL, dir.fastq.processed = NULL,
	dir.bam = NULL) {

setwd("C:/Users/assenov/DKFZ/Projects")
suppressPackageStartupMessages(library(devtools))
suppressPackageStartupMessages(load_all("echipp"))
options(echipp.verbose = TRUE)
#dataset <- echipp.initialize("D:/Datasets/ChIP-seq/2014-08-18-Finke/samples-Finke.csv", "mm10",
#	"D:/Datasets/ChIP-seq/colors.csv", TRUE)
dataset <- echipp.initialize("D:/Datasets/ChIP-seq/samples-ENCODE-Condorelli-cluster.xlsx", "mm10",
	"D:/Datasets/ChIP-seq/colors.csv", TRUE)
dir.report <- "C:/Users/assenov/DKFZ/Projects/60-Finke/reports/testme"
dir.fastq.original <- NULL
dir.fastq.processed <- "D:/Datasets/ChIP-seq/fastqc/fastq-processed"
dir.bam <- NULL

	echipp.require('RnBeads')
	logger.start(fname = NA)

	## Validate parameters
	validate.dir <- function(dir.name, param.name, optional = TRUE) {
		if (!(optional && is.null(dir.name))) {
			if (!(is.character(dir.name) && length(dir.name) == 1 && isTRUE(dir.name != ""))) {
				stop(paste("Invalid value for", param.name))
			}
			if ((optional || file.exists(dir.name)) && (!isTRUE(file.info(dir.name)[, "isdir"]))) {
				stop(paste0("Invalid value for ", param.name, "; expected ", ifelse(optional, " existing", ""),
					"directory"))
			}
		}
	}
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	validate.dir(dir.report, "dir.report", FALSE)
	validate.dir(dir.fastq.original, "dir.fastq.original")
	validate.dir(dir.fastq.processed, "dir.fastq.processed")
	validate.dir(dir.bam, "dir.bam")
	dirs.all <- c("original fastq" = dir.fastq.original, "processed fastq" = dir.fastq.processed, "BAM" = dir.bam)
	if (length(dirs.all) == 0) {
		stop("At least one of dir.fastq.original, dir.fastq.processed or dir.bam must be specified")
	}
	rm(validate.dir)

	ggplot2::theme_set(ggplot2::theme_bw())
	sample.ids <- rownames(dataset@info)
#	zips.expected <- paste0(sample.ids, "_fastqc.zip")
#	names(zips.expected) <- sample.ids

	## Initialize the report
	report <- createReport(file.path(dir.report, "combined.html"), "FastQC Results", "FastQC Results", "echipp",
		init.configuration = !file.exists(file.path(dir.report, "configuration")))
	txt <- c('This report summarizes the results of running <a ',
		'href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc">FastQC</a> on some or all of the available ',
		'<b>', length(sample.ids), '</b> samples from the study of interest. Information on how to interpret the ',
		'plots shown here is available on the dedicated <a ',
		'href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help">Help pages</a> in the FastQC web site.')
	report <- rnb.add.section(report, "Introduction", txt)

	## Parse the states of the tests
	t.levels <- c("PASS" = "passed", "WARN" = "warning", "FAIL" = "error")
	t.errors <- c("MIS" = "report is missing", "INV" = "invalid ZIP file", "MOD" = "unexpected modules",
		"RES" = "unexpected results")
	tbl.errors <- matrix(as.character(NA), nrow = length(dirs.all), ncol = length(sample.ids),
		dimnames = list(names(dirs.all), sample.ids))
	tbl.results <- list()
	for (ftype in names(dirs.all)) {
#ftype <- names(dirs.all)[1]
		tbl.results[[ftype]] <- matrix(as.character(NA), nrow = length(SUMMARY.STATES), ncol = length(sample.ids),
			dimnames = list(SUMMARY.STATES, sample.ids))
		files.present <- dir(dirs.all[ftype])
		for (sample.id in sample.ids) {
#sample.id <- sample.ids[1]
			fname <- file.path(dirs.all[ftype], paste0(sample.id, "_fastqc.zip"))
			if (file.exists(fname)) {
				fn <- paste0(sample.id, "_fastqc/summary.txt")
				states <- tryCatch(
					suppressWarnings(read.delim(unz(fname, fn), header = FALSE, stringsAsFactors = FALSE)),
					error = function(er) { NULL })
				if (is.null(states)) {
					tbl.errors[ftype, sample.id] <- "INV"
				} else if (!(ncol(states) >= 2 && setequal(states[, 2], SUMMARY.STATES))) {
					tbl.errors[ftype, sample.id] <- "MOD"
				} else if (!all(states[, 1] %in% names(t.levels))) {
					tbl.errors[ftype, sample.id] <- "RES"
				} else {
					tbl.results[[ftype]][states[, 2], sample.id] <- t.levels[states[, 1]]
				}
			} else {
				tbl.errors[ftype, sample.id] <- "MIS"
			}
		}
	}
	suppressWarnings(rm(ftype, files.present, sample.id, fname, fn, states))

	epi.plot.style <- function(i.legend) {
		ggplot2::theme(legend.position = c(1, 0.5), legend.justification = c(0, 0.5)) +
		ggplot2::theme(plot.margin = grid::unit(0.1 + c(0.1, i.legend, 0, 0), "in")) +
		ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid = ggplot2::element_blank()) +
		ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
		ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5, angle = 90))
	}

	## Display errors in processing
	i.xlab <- 2
	i.ylab <- 3
	i.legend <- 2
	if (!all(is.na(tbl.errors))) {
		txt <- c("Some of the FastQC reports were not generated or could not be processed by echipp. The following ",
			"figure shows the errors encountered while processing one or more of the FastQC reports.")
		rnb.add.paragraph(report, txt)
		dframe <- data.frame(
			x = factor(rep(sample.ids, each = nrow(tbl.errors)), levels = sample.ids),
			y = factor(rep(rownames(tbl.errors), ncol(tbl.errors)), levels = rev(rownames(tbl.errors))),
			z = factor(as.vector(tbl.errors), levels = names(t.errors)))
		dframe <- dframe[!is.na(dframe$z), ]
		levels(dframe$z) <- unname(t.errors)
		i.width <- 0.2 + i.ylab + ncol(tbl.errors) * 0.25 + i.legend
		i.height <- max(0.4 + i.xlab + (nrow(tbl.errors) + 0.0) * 0.25, 3.3)
		pp <- ggplot2::ggplot(dframe, ggplot2::aes_string(x = 'x', y = 'y', pch = 'z')) +
			ggplot2::coord_fixed() + ggplot2::geom_point() + ggplot2::labs(x = NULL, y = NULL, pch = "Error") +
			ggplot2::scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
			ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
			ggplot2::scale_shape_discrete(drop = FALSE) + epi.plot.style()
		rplot <- createReportPlot("errors", report, i.width, i.height)
		pp <- suppressWarnings(ggplot_gtable(ggplot_build(pp)))
		pp$widths[[3]] <- unit(i.ylab, "in")
		pp$heights[[length(pp$heights) - 2L]] <- unit(i.xlab, "in")
		grid.draw(pp)
		rplot <- off(rplot)
		txt <- "Failures when reading zipped FastQC reports."
		report <- rnb.add.figure(report, txt, rplot)
		rm(dframe, i.width, i.height, pp, rplot)
	}
	rm(tbl.errors)
	i.legend <- 1

	## Display test results in a heatmap
	report.plots <- list()
	col.levels <- c("#00FF00", "#FFFF00", "#800000")
	names(col.levels) <- t.levels
	mtraits <- echipp.color.mappings(dataset)
for (i in 1:length(tbl.results)) {
#i <- 1L
		tresults <- tbl.results[[i]]
		dframe <- data.frame(
			x = factor(rep(sample.ids, each = nrow(tresults)), levels = sample.ids),
			y = factor(rep(rownames(tresults), ncol(tresults)), levels = rev(rownames(tresults))),
			f = factor(as.vector(tresults), levels = unname(t.levels)))
		i.width <- 0.2 + i.ylab + ncol(tresults) * 0.25 + i.legend
		i.height <- 0.4 + i.xlab + (nrow(tresults) + 0.0) * 0.25
		for (it in 0:length(mtraits)) {
#it <- 1L
			pp <- ggplot2::ggplot(dframe, ggplot2::aes_string(x = 'x', y = 'y', fill = 'f')) +
				ggplot2::geom_tile(col = "#FFFFFF") + scale_fill_manual(na.value = "#FFFFFF", values = col.levels) +
				ggplot2::scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
				ggplot2::scale_y_discrete(drop = FALSE, expand = c(0, 0)) +
				ggplot2::labs(x = NULL, y = NULL, fill = NULL) + epi.plot.style(i.legend)
			if (it != 0) {
				s.colors <- mtraits[[it]][as.character(dataset@info[, names(mtraits)[it]])]
				for (j in 1:length(s.colors)) {
					if (!is.na(s.colors[j])) {
						pp <- pp + ggplot2::annotation_custom(
							grob = grid::rectGrob(gp = gpar(col = "#00000000", fill = s.colors[j])),
							xmin = j - 0.49, xmax = j + 0.49,
							ymin = nlevels(dframe$y) + 0.6, ymax = nlevels(dframe$y) + 0.8)
					}
				}
				rm(s.colors, j)
			}
			fname <- paste("summary", i, it, sep = "_")
			rplot <- createReportPlot(fname, report, width = i.width, height = i.height)
			pp <- suppressWarnings(ggplot2::ggplot_gtable(ggplot2::ggplot_build(pp)))
			pp$widths[[3]] <- unit(i.ylab, "in")
			pp$heights[[length(pp$heights) - 2L]] <- unit(i.xlab, "in")
			pp$layout$clip[pp$layout$name == "panel"] <- "off"
			grid.draw(pp)
			report.plots <- c(report.plots, off(rplot))
		}
	}

}
