#!/usr/bin/Rscript


# Copyright (c) 2022 University of Tartu
# Distributed under terms of the MIT Licence (see LICENCE.txt)
# Contact: Maarja Lepamets <maarja.lepamets@ut.ee>



#
# Author: Maarja Lepamets
# Date: 2019-10-03
# ============================
#
# Calculating WGS quality metric for CNVs
# detected from genotyping array. The metric is calculated
# as the fraction of CNV (in basepairs) validated by WGS data.
#
# ============================
#

# ==== Libraries
library(data.table)
library(dplyr)
library(parallel)
library(optparse)

# ==== Functions


read_samples <- function(samples) {
	if (is.null(samples)) return(NULL)
	samples <- unique(fread(samples, header=F, data.table=F)$V1)
	cat("Unique samples read from the samples file: N =", nrow(samples), "\n")
	return(samples)
}

read_pcnv <- function(pcnv) {
	pcnv <- fread(pcnv, data.table=F) 
	cat("pCNV table read...\n=== Number of unique samples:", length(unique(pcnv$Sample_Name)), "\n=== Number of pCNVs:", nrow(pcnv), "\n")
	return(pcnv)
}

read_wgscnv <- function(wgscnv) {
	wgscnv <- fread(wgscnv, data.table=F)
	wgs_probes <- wgscnv[1:3] 
	wgs_genotypes <- wgscnv[4:ncol(wgscnv)] 
	cat("WGS CNV table read... Found", ncol(wgs_probes), "header columns {", paste(names(wgs_probes), collapse = " ; "), "} and", ncol(wgs_genotypes), "samples\n")
	return(list(wgs_probes = wgs_probes, wgs_genotypes = wgs_genotypes))
}

# retain the set of mutual samples
get_filtered_samples <- function(samples, pcnv, wgscnv) {
	if (is.null(samples)) samples <- intersect(unique(pcnv$Sample_Name), names(wgscnv$wgs_genotypes))
	else samples <- intersect(samples, intersect(unique(pcnv$Sample_Name), names(wgscnv$wgs_genotypes)))
	cat("Find sample overlap: N =", length(samples), "\n")
	return(samples)
}

filter_pcnv <- function(pcnv, samples) {
	pcnv <- pcnv %>% filter(Sample_Name %in% samples)
	return(pcnv)
}

filter_wgscnv <- function(wgscnv, samples) {
	wgs_genotypes <- wgscnv$wgs_genotypes[samples]
	return(list(wgs_probes = wgscnv$wgs_probes, wgs_genotypes = wgs_genotypes))
}

# metric calculations
calculate_wgs_metric <- function (pcnv, wgscnv, cores) {

	cat("Calculate WGS metrics...\n")

	chrs <- 1:22 # only autosomes will be used
	wgscnv_probes <- wgscnv$wgs_probes
	wgscnv_genotypes <- wgscnv$wgs_genotypes

        results <- bind_rows(lapply(chrs, function(this_chr) {
		cat ("Chromosome", this_chr, "\n")

		# extract chr data
                this_chr_pcnv <- pcnv %>% filter(Chromosome == this_chr)
		if (nrow(this_chr_pcnv) == 0) return(NULL)

		this_chr_wgs_probes <- wgscnv_probes %>% filter(CHR == this_chr)
		this_chr_wgs_genotypes <- wgscnv_genotypes[wgscnv_probes$CHR == this_chr,]

		# iterate over all pCNV
                overlap <- unlist(mclapply(1:nrow(this_chr_pcnv), function(i) {

                        cnv_start <- this_chr_pcnv$Start_Position_bp[i]
                        cnv_end <- this_chr_pcnv$End_Position_bp[i]
                        cnv_sample <- this_chr_pcnv$Sample_Name[i]
			cnv_is_deletion <- this_chr_pcnv$Copy_Number[i] < 2

			# extract corresponding WGS CNV
                        cnv_probes <- which(this_chr_wgs_probes$END >= cnv_start & this_chr_wgs_probes$START <= cnv_end)
                        
			# if no overlapping WGS CNVs:
			if (length(cnv_probes) == 0) return(0)

			# otherwise extract genotypes and filter based on copy number
                        cnv_sample_genotypes <- unlist(this_chr_wgs_genotypes[cnv_probes, cnv_sample])
			if (cnv_is_deletion) cnv_probes <- cnv_probes[!is.na(cnv_sample_genotypes) & cnv_sample_genotypes < 2]
			else cnv_probes <- cnv_probes[!is.na(cnv_sample_genotypes) & cnv_sample_genotypes > 2]                        
                        
			# if now no overlapping WGS CNV of correct type:
			if (length(cnv_probes) == 0) return(0)
 
			# otherwise calculate overlap

			# 1. consider overlapping WGS CNVs only once:
                        cnv_wgs_bp <- c(this_chr_wgs_probes$START[cnv_probes], this_chr_wgs_probes$END[cnv_probes])
                        cnv_wgs_bptype <- c(rep("START", length(cnv_probes)), rep("END", length(cnv_probes)))
                        ord <- order(cnv_wgs_bp)
                        cnv_wgs_bp <- cnv_wgs_bp[ord]
                        cnv_wgs_bptype <- cnv_wgs_bptype[ord]

                        which_positions_unique <- NULL
                        tmp <- 0
                        for (j in 1:length(cnv_wgs_bptype)) {
                                tt <- cnv_wgs_bptype[j]
                                if (tmp == 0) which_positions_unique <- c(which_positions_unique, j)
                                if (tt == "START") tmp <- tmp + 1
                                if (tt == "END") tmp <- tmp - 1
                                if (tmp == 0) which_positions_unique <- c(which_positions_unique, j)
                        }
                        cnv_wgs_starts <- cnv_wgs_bp[which_positions_unique[seq(1, length(which_positions_unique), 2)]]
                        cnv_wgs_ends <- cnv_wgs_bp[which_positions_unique[seq(2, length(which_positions_unique), 2)]]

			# 2. calculate overlap between pCNV and WGS CNV:
                        overlap <- 0
                        for (j in 1:length(cnv_wgs_starts)) {
                                overlap <- overlap + (min(cnv_wgs_ends[j], cnv_end) - max(cnv_wgs_starts[j], cnv_start))
                        }
                        overlap <- overlap / (cnv_end - cnv_start + 1)
                        return(overlap)

                }, mc.cores = cores))

		this_chr_pcnv$WGS_Metric <- overlap
		this_chr_pcnv$WGS_Metric[this_chr_pcnv$Copy_Number < 2] <- -1 * this_chr_pcnv$WGS_Metric[this_chr_pcnv$Copy_Number < 2]

                return(this_chr_pcnv)

        }))
        return(results)
}


run_workflow <- function(opt) {
	# read input
	samples <- read_samples(opt$samples)
	pcnv <- read_pcnv(opt$pcnv)
	wgscnv <- read_wgscnv(opt$wgscnv)

	# filter input
	samples <- get_filtered_samples(samples, pcnv, wgscnv)
	pcnv <- filter_pcnv(pcnv, samples)
	wgscnv <- filter_wgscnv(wgscnv, samples)

	# calculate metric
	pcnv <- calculate_wgs_metric(pcnv, wgscnv, opt$cores)
	
	# write output
	write.table(pcnv, file = opt$output, row.names=F, quote=F, sep="\t")
}


# ====== Main
option_list <- list(
	make_option("--pcnv", type = "character", default = NULL,
		help = "Input table of genotyping array CNV", metavar = "character"),
	make_option("--wgscnv", type = "character", default = NULL,
		help = "Input table with WGS CNV", metavar = "character"),
	make_option(c("-s", "--samples"), type = "character", default = NULL,
		help = "File with samples to include in the metric calculations. If NULL then all overlapping samples between --pcnv and --wgscnv are used (default = NULL)", metavar = "character"),
	make_option("--cores", type = "numeric", default = 1,
		help = "Number of threads (default = 1)", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = "pcnv_wgs_metrics.tsv",
		help = "Output data table, contains genotyping array CNV and calculated WGS metrics", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$pcnv) | is.null(opt$wgscnv)) {
	print_help(opt_parser)
	stop("--pcnv and --wgscnv parameters must be provided.", call.=F)
}

run_workflow(opt)

