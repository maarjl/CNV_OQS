#!/usr/bin/Rscript


# Author: Maarja Lepamets
# Date: 2022-02-24

# Calculations of a combined omics-based metric


# ====== Libraries
library(data.table)
library(dplyr)
library(optparse)


# ===== Functions

read_metric_table <- function(input) {
	if (is.null(input)) return(NULL)

	data <- fread(input, data.table=F)
	return(data)
}

get_column_names <- function(wgs, meth, ge) {
	col_names <- NULL

	if (!is.null(wgs)) col_names <- c(col_names, names(wgs))
	if (!is.null(meth)) col_names <- c(col_names, names(meth))
	if (!is.null(ge)) col_names <- c(col_names, names(ge))

	col_names <- col_names[!(col_names %in% c("WGS_Metric", "GE_Metric", "MET_Metric"))]
	return(unique(col_names))
}

fix_empty_table <- function(table, column_names) {
	if (is.null(table)) {
		table <- data.frame(matrix(ncol = length(column_names), nrow = 0))
		colnames(table) <- column_names
	}
	return(table)
}


merge_metric_tables <- function(wgs, meth, ge) {
	data <- merge(wgs, ge, by = intersect(names(wgs), names(ge)), all = TRUE)
	data <- merge(data, meth, by = intersect(names(data), names(meth)), all = TRUE)
	return(data)
}

calculate_combined_metric <- function(full_table, num_combine) {
	metrics <- full_table[,grepl("Metric", names(full_table)), drop=F]

	combined_metric <- apply(metrics, 1, function(X) {
		w <- which.max(abs(0.5 - abs(X)))
		return(X[w])
	})
	use_metric <- rowSums(!is.na(metrics)) >= num_combine
	
	combined_metric[!use_metric] <- NA
	full_table$Combined_Metric <- combined_metric
	return(full_table)
}


run_workflow <- function(opt) {

	# read input
	wgs_tbl <- read_metric_table(opt$wgs)
	met_tbl <- read_metric_table(opt$meth)
	ge_tbl <- read_metric_table(opt$ge)

	# fix empty tables
	col_names <- get_column_names(wgs_tbl, met_tbl, ge_tbl)
	wgs_tbl <- fix_empty_table(wgs_tbl, col_names)
	met_tbl <- fix_empty_table(met_tbl, col_names)
	ge_tbl <- fix_empty_table(ge_tbl, col_names)

	# combine table
	full_tbl <- merge_metric_tables(wgs_tbl, met_tbl, ge_tbl)
	full_tbl <- calculate_combined_metric(full_tbl, opt$num_combine)

	write.table(full_tbl, file = opt$output, row.names=F, quote=F, sep="\t")
}



# ===== Main

option_list <- list(
        make_option("--wgs", type = "character", default = NULL,
                help = "CNV table with WGS Metric", metavar = "character"),
        make_option("--meth", type = "character", default = NULL,
                help = "CNV table with MET Metric", metavar = "character"),
        make_option("--ge", type = "character", default = NULL,
                help = "CNV table with GE Metric", metavar = "character"),
        make_option("--num_combine", type = "numeric", default = 1,
                help = "Number of different omics metrics required to calculate the combined metric (default = 1)", metavar = "character"),
        make_option(c("-o", "--output"), type = "character", default = "full_metrics.tsv",
                help = "Output data table, contains CNV, omics-based metrics and calculated combined metric", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$wgs) & is.null(opt$meth) & is.null(opt$ge)) {
        print_help(opt_parser)
        stop("At least one of the --wgs, --meth and --ge parameters must be provided.", call.=F)
}

run_workflow(opt)

