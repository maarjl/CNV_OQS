#!/usr/bin/Rscript

# Author: Maarja Lepamets
# Date: 2018-01-29


# This script checks which eQTLs (per gene) are not correlated with CNV.
# These will be the eQTLs that we use to correct our gene expression 
# (otherwise we correct for the CNV effect as well, which we do not want)

# ============= Libraries
library(parallel)
library(data.table)
library(dplyr)
library(tibble)
library(snpStats)
library(optparse)

# ============ Functions
# read input functions:
read_eqtl_table <- function(eqtl) {
	eqtl <- fread(eqtl, data.table=F)
	cat("Number of unique eQTLs is", length(unique(eqtl$SNP)), "for", length(unique(eqtl$Name)), "unique genes...\n")
	return(eqtl)
}

read_genes_table <- function(annotations) {
	genes <- fread(annotations, data.table=F)
	cat("Number of gene locations read: N =", nrow(genes), "\n")
	return(genes)
}

read_plink_genotypes <- function(bfile) {
	genotypes <- read.plink(bfile)$genotypes
	genotypes <- matrix(as.numeric(genotypes), ncol = ncol(genotypes), nrow = nrow(genotypes), dimnames = list(rownames(genotypes), colnames(genotypes)))

	cat("Genotypes read for", nrow(genotypes), "samples and", ncol(genotypes), "markers...\n")
	return(genotypes)
}

read_pcnv <- function(pcnv) {
        pcnv <- fread(pcnv, data.table=F)
	cat("pCNV table read...\n=== Number of unique samples:", length(unique(pcnv$Sample_Name)), "\n=== Number of pCNVs:", nrow(pcnv), "\n")
        return(pcnv)
}

# filter input functions
filter_pcnv <- function(pcnv, genotypes) {
	pcnv <- pcnv %>% filter(Sample_Name %in% rownames(genotypes))
	cat("Final number of CNV carriers: N =", length(unique(pcnv$Sample_Name)), "\n")
	return(pcnv)
}


# Functions for eQTL testing
get_overlapping_cnv_per_gene <- function(genes, pcnv) {
	overlapping_cnv_per_gene <- list()
	for (i in 1:nrow(genes)) {
		overlapping_cnv_per_gene[[genes$Name[i]]] <- pcnv %>%
			filter(Chromosome == genes$Chromosome[i], Start_Position_bp <= genes$End[i], 
				End_Position_bp >= genes$Start[i]) %>%
			select(Sample_Name, Chr = Chromosome, Start = Start_Position_bp, End = End_Position_bp, CN = Copy_Number)
	}
	return(overlapping_cnv_per_gene)
}


get_associations <- function(eqtl, genotypes, overlapping_cnv_per_gene, cores) {

	non_tagging_eqtl <- mclapply(1:nrow(eqtl), function(i) {
		snp <- eqtl$SNP[i]
		gene_name <- eqtl$Name[i]

		# cnv data
		cnv_overlapping_gene <- overlapping_cnv_per_gene[[gene_name]]
		cnv_breakpoints <- unique(c(cnv_overlapping_gene$Start, cnv_overlapping_gene$End))

		genotype_of_eqtl <- unlist(genotypes[,snp])
		if (all(genotype_of_eqtl == genotype_of_eqtl[1])) return(NULL)

		P <- NULL
		for(x in cnv_breakpoints) {
			CN <- rep(0, length(genotype_of_eqtl))
			names(CN) <- rownames(genotypes)
			CN[(cnv_overlapping_gene %>% filter(Start <= x, End >= x, CN < 2))$Sample_Name] <- -1
			CN[(cnv_overlapping_gene %>% filter(Start <= x, End >= x, CN > 2))$Sample_Name] <- 1

			min_p_value <- min(summary(lm(CN ~ genotype_of_eqtl))$coefficients[2,4], 
				summary(lm((1 * (CN > 0)) ~ genotype_of_eqtl))$coefficients[2,4], 
				summary(lm((1 * (CN < 0)) ~ genotype_of_eqtl))$coefficients[2,4], na.rm = TRUE)
			P <- c(P, min_p_value)
		}
		if (min(P) > 0.05) return(data.frame(Name = gene_name, SNP = snp))
		return(NULL)
	}, mc.cores = cores)
	non_tagging_eqtl <- do.call("rbind", non_tagging_eqtl)
	return(non_tagging_eqtl)

}


format_non_tagging_eqtl_by_gene <- function(non_tagging_eqtl) {
	unique_genes <- unique(non_tagging_eqtl$Name)

	non_tagging_eqtl_by_gene <- list()
	for(gene in unique_genes) {
		non_tagging_eqtl_by_gene[[gene]] <- as.character(non_tagging_eqtl$SNP[non_tagging_eqtl$Name == gene])
	}
	return(non_tagging_eqtl_by_gene)
}

run_workflow <- function(opt) {

	eqtl <- read_eqtl_table(opt$eqtl)
	genes <- read_genes_table(opt$gene_loc)
	genotypes <- read_plink_genotypes(opt$bfile)
	pcnv <- read_pcnv(opt$pcnv)

	pcnv <- filter_pcnv(pcnv, genotypes)

	overlapping_cnv_per_gene <- get_overlapping_cnv_per_gene(genes, pcnv)

	non_tagging_eqtl <- get_associations(eqtl, genotypes, overlapping_cnv_per_gene, opt$cores)
	non_tagging_eqtl_by_gene <- format_non_tagging_eqtl_by_gene(non_tagging_eqtl)

	saveRDS(non_tagging_eqtl_by_gene, file = opt$output)
}


# =================== Main
option_list <- list(
	make_option("--pcnv", type = "character", default = NULL,
                help = "Input table of CNV", metavar = "character"),
        make_option("--eqtl", type = "character", default = NULL,
                help = "Two-column tables of top (conditional) eQTLs", metavar = "character"),
        make_option("--gene_loc", type = "character", default = NULL,
                help = "Table with gene names and positions", metavar = "character"),
        make_option("--bfile", type = "character", default = NULL,
                help = "Prefix of genotypes of eQTL SNPs in PLINK v1 format (bed/bim/fam)", metavar = "character"),
	make_option("--cores", type = "numeric", default = 1,
                help = "Number of threads (default = 1)", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = NULL,
		help = "Output of eQTLs that do not correlate with CNV (RDS file)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$eqtl) | is.null(opt$gene_loc) | is.null(opt$bfile) | is.null(opt$pcnv) | is.null(opt$output)) {
        print_help(opt_parser)
        stop("All inputs must be provided", call.=F)
}

run_workflow(opt)

