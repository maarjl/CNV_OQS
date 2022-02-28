#!/usr/bin/Rscript

# Author: Maarja Lepamets
# Date: 2018-09-26

# Build CNV quality models

# ==== Libraries
library(plyr)
library(dplyr)
library(data.table)
library(parallel)
library(optparse)

#source("/gpfs/space/GI/GV/Projects/CNV_omics_quality_score/scripts/R-modelling_functions.R")

# ====================== Functions

# internal helping functions
.make_k_random_index_subsets <- function(N, k) {
        set.seed(1)
        ss <- sample(N)
        subsets <- list()
        for (i in 1:k) {
                subsets[[i]] <- ss[seq(i, N, k)]
        }
        return(subsets)
}

.build_interaction_terms <- function(terms1, terms2) {
        interactions <- expand.grid(terms1, terms2)
        interactions <- apply(interactions, 1, function(x) paste(x[1], x[2], sep=":"))
        return(interactions)
}

.glm_mse <- function(formula, metric, data, test_indices) {
        fit <- glm(formula, data = data[-test_indices,], family = binomial(link = "logit"))
        predictions <- predict(fit, newdata = data[test_indices,], type = "response")
        return(mean((data[test_indices,metric] - predictions)^2, na.rm = TRUE))
} 

# main modelling functions
read_metrics_table <- function(input) {
	data <- fread(input, data.table=F)
	return(data)
}

build_stepwise_logit_model <- function(data_table, metric, sample_X, cnv_X, k, k_folds, cores) {
	X_full <- c(sample_X, cnv_X)

	formula <- paste0(metric, "~")
	current_minimal_penalty <- Inf
	X_in_model <- NULL

	while (TRUE) {
		if (length(X_full) == 0) break

		# cress-validation step
		current_step_tested_models <- mclapply(X_full, function(X) {
			penalty <- NULL
			for (i in 1:k) {
				penalty <- c(penalty, .glm_mse(paste0(formula, X), metric, data_table, k_folds[[i]]))
			}
			return(data.frame(formula = paste0(formula, X), penalty = mean(penalty)))
		}, mc.cores = cores)
		current_step_tested_models <- bind_rows(current_step_tested_models)
		if (min(current_step_tested_models$penalty, na.rm = T) > current_minimal_penalty) break

		w <- which.min(current_step_tested_models$penalty)
		new_parameter <- X_full[w]
		X_full <- X_full[-w]

		# add new (interaction) parameters to X_full
		if (new_parameter %in% sample_X & any(cnv_X %in% X_in_model)) X_full <- c(X_full, .build_interaction_terms(new_parameter, X_in_model[X_in_model %in% cnv_X]))
		else if (new_parameter %in% cnv_X & any(sample_X %in% X_in_model)) X_full <- c(X_full, .build_interaction_terms(new_parameter, X_in_model[X_in_model %in% sample_X]))
		X_in_model <- c(X_in_model, new_parameter)

		formula <- paste0(current_step_tested_models$formula[w], " + ")
		current_minimal_penalty <- current_step_tested_models$penalty[w]
	}
	formula <- substr(formula, 1, nchar(formula) - 3)
	return(data.frame(formula = formula, penalty = current_minimal_penalty))
}


get_coefs_of_model <- function(formula, data) {
	fit <- glm(as.character(formula), data = data, family = binomial(link = "logit"))
	coefs <- summary(fit)$coef[,1,drop=F]
	return(coefs[,1,drop=F])
}


build_and_write_quality_model <- function(metric_table, metric_column, sample_column, cnv_column, k, cores, output) {
	if (!is.null(sample_column)) sample_column <- unlist(strsplit(sample_column, ","))
	if (!is.null(cnv_column)) cnv_column <- unlist(strsplit(cnv_column, ","))
	metric_table[,metric_column] <- abs(metric_table[,metric_column])

        deletions_table <- metric_table %>% filter(Copy_Number < 2)
        duplications_table <- metric_table %>% filter(Copy_Number > 2)

        k_folds_deletions <- .make_k_random_index_subsets(nrow(deletions_table), k)
        k_folds_duplications <- .make_k_random_index_subsets(nrow(duplications_table), k)

        del_model <- build_stepwise_logit_model(deletions_table, metric_column, sample_column, cnv_column, k, k_folds_deletions, cores)
        dup_model <- build_stepwise_logit_model(duplications_table, metric_column, sample_column, cnv_column, k, k_folds_duplications, cores)

        # calculate betas
        del_model_coefs <- get_coefs_of_model(del_model$formula, deletions_table)
	write.table(del_model_coefs, file = paste0(output, "_deletions.tsv"), col.names=F, quote=F, sep="\t")
	dup_model_coefs <- get_coefs_of_model(dup_model$formula, duplications_table)
	write.table(dup_model_coefs, file = paste0(output, "_duplications.tsv"), col.names=F, quote=F, sep="\t")

}


run_workflow <- function(opt) {
	# read input
	data_table <- read_metrics_table(opt$input)

	build_and_write_quality_model(data_table, opt$metric_column, opt$sample_column, opt$cnv_column, k=opt$k, opt$cores, opt$output)
}


# ======================= Main
option_list <- list(
        make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Full CNV table with calculated metrics, including combined metric", metavar = "character"),
	make_option("--metric_column", type = "character", default = "Combined_Metric",
		help = "Column name of the metric to predict (default: Combined_Metric)", metavar = "character"),
	make_option("--sample_column", type = "character", default = NULL,
                help = "Comma-separated list of sample-specific column names", metavar = "character"),
	make_option("--cnv_column", type = "character", default = NULL,
                help = "Comma-separated list of CNV-specific column names", metavar = "character"),
	make_option("--k", type = "numeric", default = 3,
                help = "Number of cross-validation steps (default: 3)", metavar = "character"),
	make_option("--cores", type = "numeric", default = 1,
                help = "Number of threads used for calculations (default: 1)", metavar = "character"),
        make_option(c("-o", "--output"), type = "character", default = "oqs_model",
                help = "Prefix for output models files", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (is.null(opt$input) | (is.null(opt$sample_column) & is.null(opt$cnv_column))) {
        print_help(opt_parser)
        stop("Parameters --input and at least one out of --sample_column and --cnv_column must be provided.", call.=F)
}

run_workflow(opt)
