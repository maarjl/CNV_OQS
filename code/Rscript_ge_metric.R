#!/usr/bin/Rscript


# Copyright (c) 2022 University of Tartu
# Distributed under terms of the MIT Licence (see LICENCE.txt)
# Contact: Maarja Lepamets <maarja.lepamets@ut.ee>


# Author: Maarja Lepamets
# Date: 2018-09-24


# Calculations of GE metrics



# ====== Libraries
library(dplyr)
library(tibble)
library(data.table)
library(parallel)
library(snpStats)
library(optparse)



# ====== Functions

# reading input data functions
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

read_eqtl_list <- function(eqtl) {
	if (is.null(eqtl)) return(NULL)

        eqtl <- readRDS(eqtl)
        cat("Number of unique eQTLs is", length(unique(unlist(eqtl))), "for", length(eqtl$ID), "unique genes...\n")
        return(eqtl)
}

read_genes_table <- function(annotations) {
        genes <- fread(annotations, data.table=F)
        cat("Number of gene locations read: N =", nrow(genes), "\n")
        return(genes)
}

read_plink_genotypes <- function(bfile) {
	if (is.null(bfile)) return(NULL)

        genotypes <- read.plink(bfile)$genotypes
        genotypes <- matrix(as.numeric(genotypes), ncol = ncol(genotypes), nrow = nrow(genotypes), dimnames = list(rownames(genotypes), colnames(genotypes)))

        cat("Genotypes read for", nrow(genotypes), "samples and", ncol(genotypes), "markers...\n")
        return(genotypes)
}

read_gene_expression_data <- function(ge) {
	ge_expr <- readRDS(ge) %>% as.matrix() %>% t()
	return(ge_expr)
}


# filtering functions:
get_filtered_samples <- function(samples, ge, genotypes, covariates) {
        final_samples <- rownames(ge)

	if (!is.null(genotypes)) final_samples <- intersect(final_samples, rownames(genotypes))
        if (!is.null(covariates)) final_samples <- intersect(final_samples, rownames(covariates))
        if (!is.null(samples)) final_samples <- intersect(final_samples, samples)

        cat("Final samples in gene expression data: N =", length(final_samples), "\n")
        return(final_samples)
}

filter_genes <- function(genes, dd_genes, dd_threshold) {
	if (is.null(dd_genes)) return(genes)

	gene_cnv_correlations <- fread(dd_genes, header = T, data.table=F)
        return(genes[genes$Name %in% gene_cnv_correlations$hugo_gene[gene_cnv_correlations$pearson_r >= dd_threshold],])
}

filter_pcnv <- function(pcnv, samples) {
        pcnv <- pcnv %>% filter(Sample_Name %in% samples)
        cat("Final number of CNV carriers: N =", length(unique(pcnv$Sample_Name)), "\n")
        return(pcnv)
}


filter_ge <- function(ge, samples, genes) {
        ge <- ge[samples, genes$ID]
        return(ge)
}

filter_covariates <- function(covariates, samples) {
        if (is.null(covariates)) return(NULL)
        covariates <- covariates[samples,,drop=F]
        return(covariates)
}

filter_genotypes <- function(genotypes, samples) {
	if (is.null(genotypes)) return(NULL)
	genotypes <- genotypes[samples,]
	return(genotypes)
}

# ge correction functions:
correct_for_covariates <- function(ge, covariates, cores) {
        if (is.null(covariates)) return(ge)

        residuals <- do.call("cbind", mclapply(1:ncol(ge), function(i) {
                this_res_ordered <- rep(NA, nrow(ge))
                raw_res <- lm(ge[,i] ~ ., data = covariates)$residuals
                this_res_ordered[!is.na(ge[,i]) & complete.cases(covariates)] <- raw_res
                return(this_res_ordered)
        }, mc.cores = cores))

        scaled_residuals <- apply(residuals, 2, function(x) (x - mean(x, na.rm = T)) / sd(x[!is.na(x)]))
        scaled_residuals[is.na(scaled_residuals)] <- 0

        rownames(scaled_residuals) <- rownames(ge)
        colnames(scaled_residuals) <- colnames(ge)

        return(scaled_residuals)
}


correct_for_pcs <- function(ge, npcs, cores) {
        cat("Correct for", npcs, "PCs\n")
        if (npcs == 0) return(ge)

        gram_mat <- ge %*% t(ge) / nrow(ge)
        e <- eigen(gram_mat)
        pcs <- e$vectors[,1:npcs]
        pcs <- as.data.frame(pcs)

        residuals <- do.call("cbind", mclapply(1:ncol(ge), function(i) {
                this_res_ordered <- rep(NA, nrow(ge))
                raw_res <- lm(ge[,i] ~ ., data = pcs)$residuals
                this_res_ordered[!is.na(ge[,i]) & complete.cases(pcs)] <- raw_res
                return(this_res_ordered)
        }, mc.cores = cores))

        rownames(residuals) <- rownames(ge)
        colnames(residuals) <- colnames(ge)
        return(residuals)
}


correct_for_eqtls <- function(ge, eqtls_by_gene, genotypes) {
	if (is.null(eqtls_by_gene) | is.null(genotypes)) return(ge)

	residuals <- lapply(1:ncol(ge), function(i) {
                gene_name <- colnames(ge)[i]
                if (!(gene_name %in% names(eqtls_by_gene))) return(ge[,i])

                G <- genotypes[,eqtls_by_gene[[gene_name]], drop=F] %>% as.data.frame()
                fit <- lm(ge[,i] ~ ., data = G)
                corrected_residuals <- rep(NA, nrow(G))
                corrected_residuals[complete.cases(G)] <- fit$residuals
                return(corrected_residuals)
	})
	eqtl_corrected_residuals <- do.call("cbind", residuals)
	rownames(eqtl_corrected_residuals) <- rownames(ge)
	colnames(eqtl_corrected_residuals) <- colnames(ge)
	return(eqtl_corrected_residuals)
}

# GE metric calculation functions


find_overlapping_genes_per_cnv <- function(pcnv, genes, min_overlap) {
	unique_cnvs <- pcnv %>%
		select(Chromosome, Start = Start_Position_bp, End = End_Position_bp) %>%
		unique()

	overlapping_genes_per_cnv <- list()
	for (i in 1:nrow(unique_cnvs)) {
		cnv_id <- paste0(unique_cnvs[i, ], collapse = "_")
		overlapping_genes_per_cnv[[cnv_id]] <- genes %>%			
			filter(Chromosome == unique_cnvs$Chromosome[i]) %>%
			mutate(overlap = (pmin(End, unique_cnvs$End[i]) - pmax(Start, unique_cnvs$Start[i]) + 1) / (End - Start + 1)) %>%
			filter(overlap >= min_overlap) %>%
			select(ID, Chromosome, Start, End, overlap)
	}
	
	return(overlapping_genes_per_cnv)
}

find_all_overlapping_genes <- function(overlapping_genes_per_cnv) {	
	all_overlapping_genes <- do.call("rbind", overlapping_genes_per_cnv) %>%
		select(ID, Chromosome, Start, End) %>%
		unique()
		
	return(all_overlapping_genes)
}

calculate_score_matrix <- function(pcnv, ge, all_overlapping_genes, cores) {

	score_matrix <- mclapply(1:nrow(all_overlapping_genes), function(i) {
		scores <- rep(NA, nrow(ge))
                names(scores) <- rownames(ge)

		carrier_data <- pcnv %>%
			mutate(is_near_gene = Chromosome == all_overlapping_genes$Chromosome[i] & all_overlapping_genes$Start[i] <= End_Position_bp & Start_Position_bp <= all_overlapping_genes$End[i]) %>%
			group_by(Sample_Name) %>%
			summarize(cnvs_near_gene = sum(is_near_gene)) %>%
			mutate(type = ifelse(cnvs_near_gene == 0, "noncarrier", "carrier"))
	
		carriers <- filter(carrier_data, type == "carrier")$Sample_Name
		noncarriers <- rownames(ge)[!(rownames(ge) %in% carriers)]

		expr_noncarriers <- ge[noncarriers, all_overlapping_genes$ID[i]]
		expr_carriers_std <- (ge[carriers, all_overlapping_genes$ID[i]] - mean(expr_noncarriers)) / sd(expr_noncarriers)
		cnv_type_estimate <- 1 * (expr_carriers_std > 0) - 1 * (expr_carriers_std < 0)
		scores[carriers] <- cnv_type_estimate * (2 * pnorm(abs(expr_carriers_std)) - 1)
	
		return(scores)
	}, mc.cores = cores)
	score_matrix <- do.call("cbind", score_matrix)
        colnames(score_matrix) <- all_overlapping_genes$ID
	
	return(score_matrix)
}

calculate_ge_metric <- function(pcnv, score_matrix, overlapping_genes_per_cnv) {
	pcnv$GE_Metric <- NA

	for (i in 1:nrow(pcnv)) {
		cnv_id <- paste0(pcnv[i, c("Chromosome", "Start_Position_bp", "End_Position_bp")], collapse = "_")
		overlapping_genes <- overlapping_genes_per_cnv[[cnv_id]]
	
		if (nrow(overlapping_genes) == 0)
			next
	
		scores <- score_matrix[pcnv$Sample_Name[i], overlapping_genes$ID]
	
		if (all(is.na(scores))) next
		pcnv$GE_Metric[i] <- mean(scores, na.rm = TRUE)
	}
	pcnv$GE_Metric[pcnv$Copy_Number > 2 & pcnv$GE_Metric < 0] <- 0
	pcnv$GE_Metric[pcnv$Copy_Number < 2 & pcnv$GE_Metric > 0] <- 0
	
	return(pcnv)
}


run_workflow <- function(opt) {

	# read data
	genes <- read_genes_table(opt$gene_loc)
	samples <- read_samples(opt$samples)
        pcnv <- read_pcnv(opt$pcnv)
        covariates <- read_covariates(opt$covariates)
	eqtl <- read_eqtl_list(opt$eqtl)
	genotypes <- read_plink_genotypes(opt$bfile)
	ge <- read_gene_expression_data(opt$ge)

	# filter samples
	samples <- get_filtered_samples(samples, ge, genotypes, covariates)
	pcnv <- filter_pcnv(pcnv, samples)
	covariates <- filter_covariates(covariates, samples)
	genotypes <- filter_genotypes(genotypes, samples)

	# filter genes
	genes <- filter_genes(genes, opt$dd_genes, opt$dd_cutoff)
	ge <- filter_ge(ge, samples, genes)	

	# residualisation
	ge_corrected <- correct_for_covariates(ge, covariates, opt$cores)
	ge_corrected <- correct_for_pcs(ge_corrected, opt$npcs, opt$cores)
	ge_corrected <- correct_for_eqtls(ge_corrected, eqtl, genotypes)

	# calculations
	overlapping_genes_per_cnv <- find_overlapping_genes_per_cnv(pcnv, genes, opt$gene_overlap)
	all_overlapping_genes <- find_all_overlapping_genes(overlapping_genes_per_cnv)

	score_matrix <- calculate_score_matrix(pcnv, ge_corrected, all_overlapping_genes, opt$cores)
        pcnv <- calculate_ge_metric(pcnv, score_matrix, overlapping_genes_per_cnv)

	write.table(pcnv, file = opt$output, row.names=F, quote=F, sep="\t")
}


# =================== Main

option_list <- list(
	make_option("--pcnv", type = "character", default = NULL,
                help = "Input table of CNV", metavar = "character"),
	make_option("--ge", type = "character", default = NULL,
		help = "Gene expression (Z-scores) matrix in .RDS format", metavar = "character"),
	make_option("--covariates", type = "character", default = NULL,
                help = "Gene expression covariates", metavar = "character"),
        make_option("--eqtl", type = "character", default = NULL,
                help = "eQTL list in .RDS format", metavar = "character"),
        make_option("--gene_loc", type = "character", default = NULL,
                help = "Table of gene locations", metavar = "character"),
	make_option("--samples", type = "character", default = NULL,
                help = "File with samples to include in the metric calculations. If NULL then all overlapping samples between --ge, --bfile and --covariates are used (default = NULL)", metavar = "character"),
        make_option("--bfile", type = "character", default = NULL,
                help = "Genotypes of eQTL SNPs", metavar = "character"),
	make_option("--dd_genes", type = "character", default = NULL,
                help = "Table of dosage dependent genes as found by Talevich and Hunter Shain, 2018", metavar = "character"),
	make_option("--dd_cutoff", type = "numeric", default = 0.1,
		help = "Correlation cutoff between gene expression and CNV; used to define dosage dependence (default: 0.1)", metavar = "character"),
	make_option("--gene_overlap", type = "numeric", default = 0.8,
		help = "Minimum required overlap between gene and CNV (default: 0.8)", metavar = "character"),
	make_option("--cores", type = "numeric", default = 1,
                help = "Number of threads (default = 1)", metavar = "character"),
        make_option("--npcs", type = "numeric", default = 0,
                help = "Number of PCs using in the residualisation step (default: 0)", metavar = "character"),
        make_option(c("-o", "--output"), type = "character", default = "ge_metric.tsv",
                help = "Output directory of gene expression metrics", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$pcnv) | is.null(opt$ge) | is.null(opt$gene_loc)) {
        print_help(opt_parser)
        stop("--pcnv, --ge and --gene_loc parameters must be provided.", call.=F)
}

run_workflow(opt)


