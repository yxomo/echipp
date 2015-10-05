########################################################################################################################
## overview.R
## created: 2014-12-04
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions generating overview of a dataset.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

echipp.report.low.dimension <- function(dataset, dir.report) {
	#setwd("C:/Users/assenov/DKFZ/Projects")
	#suppressPackageStartupMessages(library(echipp))
	#options(echipp.verbose = TRUE)
	#dataset <- echipp.initialize("D:/Datasets/ChIP-seq/samples-ENCODE-Condorelli-Windows.csv", "mm10",
	#	"D:/Datasets/ChIP-seq/colors.csv", TRUE)
	#dir.report <- "C:/Users/assenov/DKFZ/Projects/60-Finke/reports/testme"

	echipp.require('RnBeads')
	logger.start(fname = NA)
	## Validate parameters
	if (!inherits(dataset, "EchippSet")) {
		stop("Invalid value for dataset")
	}
	echipp.validate.dir(dir.report, "dir.report", FALSE)

	## Initialize the report
	ggplot2::theme_set(ggplot2::theme_bw())
	sample.ids <- rownames(dataset@info)

	## Initialize the report
	report <- createReport(file.path(dir.report, "exploratory.html"), "Overview", "Overview", "echipp",
		init.configuration = !file.exists(file.path(dir.report, "configuration")))
	txt <- c('This report presents an overview of the <b>', length(sample.ids), '</b> samples of the study of ',
		'interest.')
	report <- rnb.add.section(report, "Introduction", txt)

	chrom.lengths <- sapply(echipp:::echipp.load.assembly(dataset@genome)$Bands, function(x) { max(x[, 1]) })
	gtiles <- GenomicRanges::tileGenome(chrom.lengths, tilewidth = 1000, cut.last.tile.in.chrom = TRUE)

	bam.params <- ScanBamParam(
		flag = scanBamFlag(isSecondaryAlignment = FALSE),
		what = c("qname", "rname", "strand", "pos", "mapq"),
		which = gtiles)
	fnames <- sapply(sample.ids, function(id) { echipp.get.files(dataset, id, "bam") })


	fname <- fnames[1]
	fname.index <- sub("\\.bam$", ".bai", fname)
	if (!file.exists(fname.index)) {
		indexBam(fname)
	}
	result <- countBam(fname, param = bam.params)
	files.summary <- do.call(rbind, lapply(fnames, function(x) { countBam(x, param = bam.params) }))
	
}
