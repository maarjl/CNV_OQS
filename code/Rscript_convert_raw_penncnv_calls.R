#!/usr/bin/Rscript


# Author: Maarja Lepamets
# Date: 2022-03-01


# Preparations of PennCNV output for CNV quality estimations

# ==== Libraries
library(data.table)
library(dplyr)
library(stringr)
library(parallel)
library(optparse)

options(scipen=999)

# === Functions

read_cnv_data <- function(cnv) {
	data <- fread(cnv, data.table=F, header=F, stringsAsFactors=F, sep=NULL)$V1
	return(data)
}

read_sample_data <- function(sample) {
	data <- fread(sample, data.table=F, stringsAsFactors=F, sep="\t")
	data$File <- substr(data$File, max(str_locate_all(data$File, "[./]")[[1]][,"start"]) + 1, nchar(data$File))
	return(data)
}

parse_cnv_data <- function(data, cores) {
	line_pattern <- "^chr([X|\\d]{1,2}):(\\d+)-(\\d+)\\s*numsnp=(\\d+)\\s*length=([,\\d]+)\\s*state\\d,cn=(\\d)\\s*([./\\w]+).*conf=[-]?([.\\d]+)"
	
	parsed_data <- do.call("rbind", mclapply(data, function(line) {
		capture_groups <- str_match(line, line_pattern)[,-1]
		# fix sample name
		capture_groups[7] <- substr(capture_groups[7], max(str_locate_all(capture_groups[7], "[./]")[[1]][,"start"]) + 1, nchar(capture_groups[7]))
		# fix length
		capture_groups[5] <- gsub(",", "", capture_groups[5])
		return(capture_groups)
	}, mc.cores = cores))

	parsed_data <- as.data.frame(parsed_data, stringsAsFactors=F) 
	names(parsed_data) <- c("Chromosome", "Start_Position_bp", "End_Position_bp", "No_Probes", "Length_bp", "Copy_Number", "Sample_Name", "Max_Log_BF")

	return(parsed_data)
}

merge_data <- function(cnv_data, sample_data) {
	full_data <- merge(cnv_data, sample_data, by.x="Sample_Name", by.y="File", all.x=T, sort=F)
	return(full_data)
}

add_additional_column <- function(full_data, arraysize) {
	full_data$Length_per_Probe <- as.numeric(full_data$Length_bp) / as.numeric(full_data$No_Probes)
	if (!is.null(arraysize)) full_data$NumCNV_bin <- 1 * (full_data$NumCNV / arraysize > 3.8e-5)
	full_data$WF <- abs(full_data$WF)
	return(full_data)
}

filter_data <- function(full_data, numcnv_thres, max_length_thres) {
	if (!is.null(numcnv_thres)) {
		full_data <- full_data %>% filter(NumCNV <= numcnv_thres)
	}
	if (!is.null(max_length_thres)) {
		max_len_tbl <- full_data %>% group_by(Sample_Name) %>%
			summarise(M = max(as.numeric(Length_bp))) %>%
			filter(M > max_length_thres)
		full_data <- full_data %>% filter(!(Sample_Name %in% max_len_tbl$Sample_Name))
	}
	return(full_data)
}

run_workflow <- function(opt) {
	cnv_data <- read_cnv_data(opt$cnv)
	cnv_data <- parse_cnv_data(cnv_data, opt$cores)
	sample_data <- read_sample_data(opt$samples)

	full_data <- merge_data(cnv_data, sample_data)
	full_data <- add_additional_column(full_data, opt$arraysize)
	full_data <- filter_data(full_data, opt$max_numcnv, opt$max_len)

	write.table(full_data, file = opt$output, row.names=F, quote=F, sep="\t")
}

# === Main

option_list <- list(
        make_option("--cnv", type = "character", default = NULL,
                help = "Cleaned PennCNV output (output from clean_cnv.pl script)", metavar = "character"),
        make_option("--samples", type = "character", default = NULL,
                help = "Sample-specific quality metrics from PennCNV (output from filter_cnv.pl, -qcsumout flag)", metavar = "character"),
        make_option("--cores", type = "numeric", default = 1,
                help = "Number of threads used in formatting", metavar = "character"),
        make_option("--max_numcnv", type = "numeric", default = NULL,
                help = "If given, all samples with number of CNVs larger will be excluded (default = NULL)", metavar = "character"),
	make_option("--max_len", type = "numeric", default = NULL,
                help = "If given, all samples with longer CNVs will be excluded (default = NULL)", metavar = "character"),
	make_option("--arraysize", type = "numeric", default = NULL,
                help = "If given, additional array size-specific parameters will be calculated (default = NULL)", metavar = "character"),
        make_option(c("-o", "--output"), type = "character", default = "pcnv_table.tsv",
                help = "Output CNV table", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$cnv) | is.null(opt$samples)) {
        print_help(opt_parser)
        stop("Parameters --cnv and --samples must be provided.", call.=F)
}

run_workflow(opt)



