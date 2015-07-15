########################################################################################################################
## bands.R
## created: 2015-01-27
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## 
########################################################################################################################

## F U N C T I O N S ###################################################################################################

echipp.half.ellipse.coords <- function(xmin, xmax, ymin, ymax, type = 1L, n = 60L) {
	xc <- (xmin + xmax) / 2
	yc <- ifelse(type == 1, ymin, ymax)
	xv <- seq(xmin, xmax, length.out = n)
	yv <- sqrt((ymax - ymin)^2 * (1 - (4 * (xv - xc)^2 / (xmax - xmin)^2)))
	if (type == 1) {
		yv <- yv + ymin
	} else {
		xv <- rev(xv)
		yv <- rev(ymax - yv)
	}
 	data.frame(x = xv, y = yv)
}

########################################################################################################################

echipp.fill.up.dframe <- function(dframe) {
	dframe$xmin <- dframe$xmax <- dframe$ymin <- dframe$ymax <- 1L
	if (!("stain" %in% colnames(dframe))) {
		dframe$stain <- "neg"
	}
	dframe
}

########################################################################################################################

#' echipp.chrom.plot
#' 
#' Creates a karyotype plot of the specified chromosome.
#' 
#' @param assembly.code Genome assembly of interest specified as a \code{character}.
#' @param chrom         Chromosome of interest.
#' @param focus.region  Two-element vector of non-negative \code{integer}s storing the genomic region to show a red
#'                      rectangle around. Setting this parameter to \code{NULL} (default) disables focusing on a region.
#' @return The initialized plot as a \code{ggplot} instance.
#' 
#' @seealso \code{\link{echipp.assemblies}} for a list of all supported genome assemblies.
#' 
#' @author Yassen Assenov
#' @export
echipp.chrom.plot <- function(assembly.code, chrom, focus.region = NULL) {

	## Validate parameters
	if (!(is.character(assembly.code) && length(assembly.code) == 1 && isTRUE(assembly.code != ""))) {
		stop("invalid value for assembly.code")
	}
	if (!(is.character(chrom) && length(chrom) == 1 && isTRUE(chrom != ""))) {
		stop("invalid value for chrom")
	}
	if (!is.null(focus.region)) {
		if (is.double(focus.region) && isTRUE(all(focus.region == as.integer(focus.region)))) {
			focus.region <- as.integer(focus.region)
		}
		if (!(is.integer(focus.region) && length(focus.region) == 2 && isTRUE(all(focus.region >= 0L)))) {
			stop("invalid value for focus.region")
		}
		focus.region <- sort(as.vector(focus.region))
	}

	## Load the chromosomal band colors and genome assembly information
	fname <- system.file("data/band.colors.RDS", package = "echipp")
	band.colors <- tryCatch(readRDS(fname), error = function(er) { NULL })
	if (is.null(band.colors)) {
		stop(paste("INTERNAL ERROR: Could not load", fname))
	}
	genome.annotation <- echipp.load.assembly(assembly.code)
	if (!(chrom %in% names(genome.annotation$Bands))) {
		stop("unsupported chromosome")
	}
	tbl <- genome.annotation$Bands[[chrom]]
	rm(genome.annotation)

	## Construct all geometric figures for the plot
	N <- nrow(tbl)
	in.bands <- rep(FALSE, N)
	dfr.border <- NULL
	dfr.border.rev <- NULL
	dfr.polygons <- list()
	for (i in 1:(N - 1)) {
		if (i == 1L) {
			dfr <- echipp.half.ellipse.coords(0, 1, tbl[1, 1], tbl[2, 1], 2L)
			dfr.border <- rbind(dfr.border, dfr)
			dfr.border.rev <- dfr[1, ]
			if (tbl[i, "stain"] != "neg") {
				dfr$stain <- tbl[i, "stain"]
				dfr.polygons[[length(dfr.polygons) + 1]] <- echipp.fill.up.dframe(dfr)
			}
		} else {
			x.last <- ifelse(tbl[i - 1, "stain"] == "acen", 0.2, 0)
			if (i == N - 1) {
				dfr <- echipp.half.ellipse.coords(0, 1, tbl[i, 1], tbl[N, 1], 1L)
				dfr.border <- rbind(dfr.border, dfr)
				dfr.border <- rbind(dfr.border, dfr.border.rev[nrow(dfr.border.rev):1, ])
				if (tbl[i, "stain"] != "neg") {
					dfr$stain <- tbl[i, "stain"]
					dfr.polygons[[length(dfr.polygons) + 1]] <- echipp.fill.up.dframe(dfr)
				}
			} else if (tbl[i, "stain"] == "acen") {
				x.next <- ifelse(isTRUE(tbl[i + 1, "stain"] == "acen"), 0.2, 0)
				xs <- c(x.last, x.next, 1 - x.next, 1 - x.last)
				ys <- c(tbl[i, 1], tbl[i + 1, 1], tbl[i + 1, 1], tbl[i, 1])
				dfr <- data.frame(x = xs, y = ys, stain = tbl[i, "stain"])
				dfr.polygons[[length(dfr.polygons) + 1]] <- echipp.fill.up.dframe(dfr)
				dfr <- data.frame(x = xs[1:2], y = ys[1:2])
				dfr.border <- rbind(dfr.border, dfr)
				dfr <- data.frame(x = xs[4:3], y = ys[4:3])
				dfr.border.rev <- rbind(dfr.border.rev, dfr)
				rm(x.next, xs, ys)
			} else {
				in.bands[i] <- TRUE
			}
		}
	}
	rm(dfr.border.rev, i, dfr)

	## Create the karyotype plot
	coords <- range(tbl[, 1])
	dframe <- data.frame(
		x = 0,
		y = 1L,
		xmin = 0,
		xmax = 1,
		ymin = tbl[in.bands, 1],
		ymax = tbl[which(in.bands) + 1L, 1] - 1L,
		stain = tbl[in.bands, "stain"])
	pp <- ggplot2::ggplot(dframe) +
		ggplot2::aes(x = x, y = y, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = stain) +
		ggplot2::geom_rect()
	for (dfr in dfr.polygons) {
		pp <- pp + ggplot2::geom_polygon(data = dfr)
	}
	pp <- pp + ggplot2::geom_path(data = echipp.fill.up.dframe(dfr.border))
	if (!is.null(focus.region)) {
		## Add rectangle and lines focusing on a region
		dfr <- echipp.fill.up.dframe(data.frame(x = c(0, 1, 1, 0), y = rep(focus.region, each = 2)))
		dfr <- rbind(dfr, dfr[1, ])
		pp <- pp + ggplot2::geom_path(data = dfr, lwd = 2, color = "red")
		dfr <- data.frame(x = c(1, 2, NA, 1, 2), y = c(focus.region[1], coords[1], NA, focus.region[2], coords[2]))
		dfr <- echipp.fill.up.dframe(dfr)
		pp <- pp + ggplot2::geom_path(data = dfr, color = "red", lty = 3L)
	}
	pp <- pp + ggplot2::scale_x_continuous(limits = c(0, 2 - as.integer(is.null(focus.region))), expand = c(0, 0)) +
		ggplot2::scale_y_continuous(limits = coords, expand = c(0, 0)) +
		ggplot2::scale_fill_manual(values = band.colors) +
		ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
		ggplot2::theme(panel.background = ggplot2::element_blank(), panel.border = ggplot2::element_blank()) +
		ggplot2::theme(panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_blank()) +
		ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank()) +
		ggplot2::theme(axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank()) +
		ggplot2::theme(plot.margin = grid::unit(0.1 + c(0, 0, 0, 0), "in"), legend.position = "none")
	return(pp)
}
