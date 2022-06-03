#!/usr/bin/Rscript


# Copyright (c) 2022 University of Tartu
# Distributed under terms of the MIT Licence (see LICENCE.txt)
# Contact: Maarja Lepamets <maarja.lepamets@ut.ee>


# Author: Maarja Lepamets
# Date: 2022-02-14


# Calculations of methylation-based CNV
# quality metrics.


# ======================= Libraries
library(data.table)
library(dplyr)
library(parallel)
library(optparse)

# ======================= Functions

# Functions for reading data
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

read_covariates <- function(covariates) {
	if (is.null(covariates)) return (NULL)
	covariates <- fread(covariates, data.table=F)

	samples <- covariates$Sample_Name
	covariates <- covariates[,-1,drop=F]
	rownames(covariates) <- samples

	cat("Number of samples in covariates file: N =", nrow(covariates), "\n")
	return(covariates)
}

read_and_filter_methylation_data <- function(meth) {
	meth_data <- readRDS(meth)
	meth_matrix <- meth_data$M + meth_data$U
	meth_p <- meth_data$detection_p[rownames(meth_matrix), colnames(meth_matrix)]
	meth_matrix[meth_data$detection_p > 1e-16] <- NA
	meth_matrix <- scale(meth_matrix) %>% t()
	meth_matrix[is.na(meth_matrix)] <- 0
	cat("Methylation data read...\n")
	return(meth_matrix)
}

read_cpg_probes <- function(cpg) {
	cpg <- fread(cpg, data.table=F)
	cat("Number of CpG sites: N =", nrow(cpg), "\n")
	return(cpg)
}



# Functions for filtering data
get_filtered_samples <- function(samples, meth, covariates) {
	final_samples <- rownames(meth)
	
	if (!is.null(covariates)) final_samples <- intersect(final_samples, rownames(covariates))
	if (!is.null(samples)) final_samples <- intersect(final_samples, samples)

	cat("Final samples in methylation data: N =", length(final_samples), "\n")
	return(final_samples)
}


filter_pcnv <- function(pcnv, samples) {
        pcnv <- pcnv %>% filter(Sample_Name %in% samples)
	cat("Final number of CNV carriers: N =", length(unique(pcnv$Sample_Name)), "\n")
        return(pcnv)
}


filter_meth <- function(meth, samples) {
	meth <- meth[samples,]
	return(meth)
}

filter_covariates <- function(covariates, samples) {
	if (is.null(covariates)) return(NULL)
	covariates <- covariates[samples,,drop=F]
	return(covariates)
}


# Functions for residualisation
correct_for_covariates <- function(meth, covariates, cores) {
	if (is.null(covariates)) return(meth)

	residuals <- do.call("cbind", mclapply(1:ncol(meth), function(i) {
		this_res_ordered <- rep(NA, nrow(meth))
		raw_res <- lm(meth[,i] ~ ., data = covariates)$residuals
		this_res_ordered[!is.na(meth[,i]) & complete.cases(covariates)] <- raw_res
		return(this_res_ordered)
	}, mc.cores = cores))

	scaled_residuals <- apply(residuals, 2, function(x) (x - mean(x, na.rm = T)) / sd(x[!is.na(x)]))
	scaled_residuals[is.na(scaled_residuals)] <- 0

	rownames(scaled_residuals) <- rownames(meth)
	colnames(scaled_residuals) <- colnames(meth)

	return(scaled_residuals)
}


correct_for_pcs <- function(meth, npcs, cores) {
	cat("Correct for", npcs, "PCs\n")
	if (npcs == 0) return(meth)

	gram_mat <- meth %*% t(meth) / nrow(meth)
	e <- eigen(gram_mat)
	pcs <- e$vectors[,1:npcs]
	pcs <- as.data.frame(pcs)

	residuals <- do.call("cbind", mclapply(1:ncol(meth), function(i) {
                this_res_ordered <- rep(NA, nrow(meth))
                raw_res <- lm(meth[,i] ~ ., data = pcs)$residuals
                this_res_ordered[!is.na(meth[,i]) & complete.cases(pcs)] <- raw_res
                return(this_res_ordered)
        }, mc.cores = cores))

	rownames(residuals) <- rownames(meth)
	colnames(residuals) <- colnames(meth)
	return(residuals)
}


find_overlapping_probes_per_cnv <- function(pcnv, cpg, cores) {
        unique_cnvs <- pcnv %>%
                select(Chromosome, Start = Start_Position_bp, End = End_Position_bp) %>%
                unique()

        overlapping_probes_per_cnv <- mclapply(1:nrow(unique_cnvs), function(i) {
                return(cpg %>% filter(Chromosome == unique_cnvs$Chromosome[i], Position >= unique_cnvs$Start[i], Position <= unique_cnvs$End[i]))
        }, mc.cores = cores)

        names(overlapping_probes_per_cnv) <- apply(unique_cnvs, 1, paste0, collapse = "_")
        return(overlapping_probes_per_cnv)
}


find_all_overlapping_probes <- function(overlapping_probes_per_cnv) {
        all_overlapping_probes <- do.call("rbind", overlapping_probes_per_cnv) %>%
                unique()
        return(all_overlapping_probes)
}


calculate_score_matrix <- function(pcnv, meth, all_overlapping_probes, cores) {

        score_matrix <- mclapply(1:nrow(all_overlapping_probes), function(i) {
                scores <- rep(NA, nrow(meth))
                names(scores) <- rownames(meth)

                carrier_data <- pcnv %>%
                        mutate(is_on_probe = Chromosome == all_overlapping_probes$Chromosome[i] & Start_Position_bp <= all_overlapping_probes$Position[i] & End_Position_bp >= all_overlapping_probes$Position[i]) %>%
                        group_by(Sample_Name) %>%
                        summarize(cnvs_on_probe = sum(is_on_probe)) %>%
                        mutate(type = ifelse(cnvs_on_probe == 0, "noncarrier", "carrier"))

                carriers <- filter(carrier_data, type == "carrier")$Sample_Name
                #noncarriers <- filter(carrier_data, type == "noncarrier")$Sample_Name
		noncarriers <- rownames(meth)[!(rownames(meth) %in% carriers)]

                meth_noncarriers <- meth[noncarriers, all_overlapping_probes$ID[i]]
                meth_carriers_std <- (meth[carriers, all_overlapping_probes$ID[i]] - mean(meth_noncarriers, na.rm = T)) / sd(meth_noncarriers[!is.na(meth_noncarriers)])
                cnv_type_estimate <- 1 * (meth_carriers_std > 0) - 1 * (meth_carriers_std < 0)
                scores[carriers] <- cnv_type_estimate * (2 * pnorm(abs(meth_carriers_std)) - 1)

                return(scores)
        }, mc.cores = cores)
        score_matrix <- do.call("cbind", score_matrix)
        colnames(score_matrix) <- all_overlapping_probes$ID
        return(score_matrix)
}

calculate_met_metric <- function(pcnv, score_matrix, overlapping_probes_per_cnv) {
        pcnv$MET_Metric <- NA
        for (i in 1:nrow(pcnv)) {
                cnv_id <- paste0(pcnv[i, c("Chromosome", "Start_Position_bp", "End_Position_bp")], collapse = "_")
                overlapping_probes <- overlapping_probes_per_cnv[[cnv_id]]

                if (nrow(overlapping_probes) == 0)
                        next

                scores <- score_matrix[pcnv$Sample_Name[i], overlapping_probes$ID]

                if(all(is.na(scores))) next
                pcnv$MET_Metric[i] <- mean(scores, na.rm = TRUE)

        }
        pcnv$MET_Metric[pcnv$Copy_Number > 2 & pcnv$MET_Metric < 0] <- 0
        pcnv$MET_Metric[pcnv$Copy_Number < 2 & pcnv$MET_Metric > 0] <- 0

        return(pcnv)
}


run_workflow <- function(opt) {
	# read input
	cpg <- read_cpg_probes(opt$cpg)
	samples <- read_samples(opt$samples)
	pcnv <- read_pcnv(opt$pcnv)
	covariates <- read_covariates(opt$covariates)
	meth <- read_and_filter_methylation_data(opt$meth)

	# filter
	samples <- get_filtered_samples(samples, meth, covariates)
	pcnv <- filter_pcnv(pcnv, samples)
	meth <- filter_meth(meth, samples)
	covariates <- filter_covariates(covariates, samples)

	# Residual calculations
	meth_corrected <- correct_for_covariates(meth, covariates, opt$cores)
	meth_corrected <- correct_for_pcs(meth_corrected, opt$npcs, opt$cores)

	# find overlapping CpG probes per CNV
	overlapping_probes_per_cnv <- find_overlapping_probes_per_cnv(pcnv, cpg, opt$cores)
	all_overlapping_probes <- find_all_overlapping_probes(overlapping_probes_per_cnv)

	# Metric calculations
	score_matrix <- calculate_score_matrix(pcnv, meth_corrected, all_overlapping_probes, opt$cores)
	pcnv <- calculate_met_metric(pcnv, score_matrix, overlapping_probes_per_cnv)

	# Write output
	write.table(pcnv, file = opt$output, row.names=F, quote=F, sep="\t")
}


# ========================== Main

option_list <- list(
        make_option("--pcnv", type = "character", default = NULL,
                help = "Input table of CNV", metavar = "character"),
	make_option("--meth", type = "character", default = NULL,
		help = "Input file (.RDS format) containing methylation intensities", metavar = "character"),
	make_option("--covariates", type = "character", default = NULL,
		help = "Covariates table (first column: Sample_Name, followed by covariate columns", metavar = "character"),
	make_option("--samples", type = "character", default = NULL,
		help = "File with samples to include in the metric calculations. If NULL then all overlapping samples between --meth and --covariates are used (default = NULL)", metavar = "character"),
	make_option("--cpg", type = "character", default = NULL,
		help = "Input table of CpG site locations (column names: Chromosome, Position)", metavar = "character"),
	make_option("--cores", type = "numeric", default = 1,
		help = "Number of threads (default = 1)", metavar = "character"),
	make_option("--npcs", type = "numeric", default = 0,
		help = "Number of PCs using in the residualisation step (default: 0)", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = NULL,
		help = "Output data table, contains CNV and calculated MET metrics", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$pcnv) | is.null(opt$meth) | is.null(opt$cpg)) {
        print_help(opt_parser)
        stop("--pcnv, --meth and --cpg parameters must be provided.", call.=F)
}

run_workflow(opt)


