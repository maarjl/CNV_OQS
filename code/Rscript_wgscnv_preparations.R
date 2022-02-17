#!/usr/bin/Rscript


#
# Author: Maarja Lepamets
# Date: 2022-02-09
# ============================
#
# Preprocessing and filtering WGS CNV data
#
# ============================
#

# ====== Libraries
library(data.table)
library(stringr)
library(dplyr)
library(optparse)

# ===== Functions

read_data <- function(input) {
	data <- fread(input, data.table=F)

	# first row for sample names
	first_row <- unlist(data[1,4:ncol(data)])
	samples <- str_extract(first_row, "[^:]+")

	names(data) <- c(c("CHR", "START", "END", samples))
	return(data)
}


.recognise_format_order <- function(cell) {
	tags <- c("CN=", "CNQ=", "FT=")

	cn_loc <- str_locate(cell, "CN=")
	cnq_loc <- str_locate(cell, "CNQ=")
	ft_loc <- str_locate(cell, "FT=")
	locs <- rbind(cn_loc, cnq_loc, ft_loc)
	o <- order(locs[,1], na.last = NA)
	tags <- tags[o]

	return(tags)
}

.extract_tag_data <- function(data_matrix, tag_pattern, tags_list, as.numeric = F) {
	end_pattern <- ":"
	if (tag_pattern == tags_list[length(tags_list)]) end_pattern <- "$"

	data_values <- str_match(data_matrix, paste0(tag_pattern, "(.*?)", end_pattern))[,2]
	if (as.numeric) data_values <- as.numeric(data_values)
	tag_data <- matrix(data_values, ncol = ncol(data_matrix), nrow = nrow(data_matrix), byrow = FALSE, dimnames = list(NULL, colnames(data_matrix)))
	return(tag_data)
}

format_data <- function(data, cnq_threshold) {
	data_cells <- as.matrix(data[,4:ncol(data)])
	tags <- .recognise_format_order(data_cells[1,1])

	data_cn <- .extract_tag_data(data_cells, "CN=", tags)
	if (is.null(cnq_threshold)) {
		cat("No CNQ threshold provided, using FT tag for quality filtering.\n")
		if (!("FT=" %in% tags)) {
			cat("Error: no FT tag provided.\n")
			quit(status = 100)
		}
		data_ft <- .extract_tag_data(data_cells, "FT=", tags)
		data_filter <- data_ft == "LQ"
	} else {
		cat("Using CNQ tag and threshold for quality filtering.\n")
		data_cnq <- .extract_tag_data(data_cells, "CNQ=", tags, as.numeric = T)
		data_filter <- data_cnq < cnq_threshold
	}
	data_cn[data_filter] <- NA

	return(data_cn)
}

get_filter_indices <- function(data, missing_threshold) {
	mis_per_site <- apply(data, 1, function(X) return(sum(is.na(X))/ncol(data)))
	idx <- which(mis_per_site < missing_threshold)
	return(idx)
}


run_workflow <- function(opt) {
	data <- read_data(opt$input)

	# extract columns with site positions
	data_pos <- data[1:3]
	# extract and format copy number genotypes
	data_cn <- format_data(data, opt$cnq)

	# calculate and filter based on missingness
	filt_idx <- get_filter_indices(data_cn, opt$missing)

	# format final data table
	data_pos <- data_pos[filt_idx,]
	data_cn <- data_cn[filt_idx,]
	final_data <- cbind(data_pos, data_cn)
	cat("Number of CNV regions in final output =", nrow(final_data), ".\n")

	write.table(final_data, file = opt$output, row.names=F, quote=F, sep="\t")
}


# ===== Main


option_list <- list(
	make_option(c("-i", "--input"), type = "character", default = NULL,
		help = "Input table of WGS copy number data converted from VCF", metavar = "character"),
	make_option("--cnq", type = "numeric", default = NULL,
		help = "Copy number quality threshold, under which the CN values will be set to missing. If 'NULL' FT tag is used for filtering.", metavar = "character"),
	make_option("--missing", type = "numeric", default = 0.1,
		help = "Threshold value for missingness, over which the CNV regions will be eliminated.", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = "wgs_cnv_formatted.tsv",
		help = "Name of the output file.", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
	print_help(opt_parser)
	stop("--input must be provided.", call.=F)
}

run_workflow(opt)
