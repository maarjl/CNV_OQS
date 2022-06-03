#!/usr/bin/Rscript


# Copyright (c) 2022 University of Tartu
# Distributed under terms of the MIT Licence (see LICENCE.txt)
# Contact: Maarja Lepamets <maarja.lepamets@ut.ee>


# Author: Maarja Lepamets
# Date: 2021-05-04

# Predictions of CNV quality scores


# ====================== Libraries
library(data.table)
library(dplyr)
library(optparse)

# ==================== Functions
read_pcnv_data <- function(pcnv) {
	pcnv <- fread(pcnv, data.table=F)
	return(pcnv)
}

read_omics_model <- function(model_file) {
	model <- read.table(model_file, header = F)
	coef <- cbind(model[,2], rep(0, nrow(model)))
	rownames(coef) <- as.character(model[,1])
	colnames(coef)    <- c('Estimate', 'Pr(>|z|)')
	return(coef)
}


score_cnv <- function(pcnv, model) {
	score <- rep(model['(Intercept)', 'Estimate'], nrow(pcnv))
	for (par_idx in 2:nrow(model)) {
		par_name <- rownames(model)[par_idx]

		if (grepl(":", par_name)) {
			par_names <- unlist(strsplit(par_name, ":"))
			score <- score + model[par_name, 'Estimate'] * pcnv[,par_names[1]] * pcnv[,par_names[2]]
		} else {
			score <- score + model[par_name, 'Estimate'] * pcnv[,par_name]
		}
	}
	score <- 1 / (1 + exp(-score))
	return(score)
}

add_omics_score <- function(pcnv, del_model, dup_model) {

	pcnv$Omics_Score <- rep(0, nrow(pcnv))
	pcnv$Omics_Score[pcnv$Copy_Number < 2] <- -1 * score_cnv(pcnv %>% filter(Copy_Number < 2), del_model)
	pcnv$Omics_Score[pcnv$Copy_Number > 2] <- score_cnv(pcnv %>% filter(Copy_Number > 2), dup_model)

	return(pcnv)
}


run_workflow <- function(opt) {

	# == Data
	pcnv <- read_pcnv_data(opt$pcnv)
	deletions_model <- read_omics_model(opt$del_model)
	duplications_model <- read_omics_model(opt$dup_model)

	pcnv <- add_omics_score(pcnv, deletions_model, duplications_model)
	write.table(pcnv, file = opt$output, row.names=F, quote=F, sep="\t")
}

# ========================== Main

option_list <- list(
	make_option("--pcnv", type = "character", default = NULL,
		help = "CNV input table containing predictor variables in columns", metavar = "character"),
	make_option("--del_model", type = "character", default = NULL,
		help = "Model file for deletion quality calculations", metavar = "character"),
	make_option("--dup_model", type = "character", default = NULL,
		help = "Model file for duplication quality calculations", metavar = "character"),
        make_option(c("-o", "--output"), type = "character", default = "pcnv_oqs_predictions.tsv",
                help = "Output CNV table with predicted quality", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$pcnv) | is.null(opt$del_model) | is.null(opt$dup_model)) {
	print_help(opt_parser)
        stop("Parameters --pcnv, --del_model and --dup_model are required.", call.=F)
}

run_workflow(opt)




