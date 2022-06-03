#!/usr/bin/Rscript

# Copyright (c) 2022 University of Tartu
# Distributed under terms of the MIT Licence (see LICENCE.txt)
# Contact: Maarja Lepamets <maarja.lepamets@ut.ee>



# Author: Maarja Lepamets
# Date: 2022-02-14


# Script for methylation data setup from .idat files

# ================ Libraries
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(optparse)

# ================== Functions
get_RGset <- function(samples, bg_correction = FALSE) {	
	
	RGset <- read.metharray(as.character(samples$Path), verbose = TRUE)
	if (bg_correction)
		RGset <- bgcorrect.illumina(RGset)  # Illumina background subtraction
	
	return(RGset)
}

get_signals_type <- function(RGset, probe_type = c("II", "I-Green", "I-Red"), M_colour, U_colour, sample_names) {
	probe_type <- match.arg(probe_type)
	probes <- getProbeInfo(RGset, type = probe_type)
	
	M <- M_colour[probes[, ifelse(probe_type == "II", "AddressA", "AddressB")], , drop = FALSE]
	rownames(M) <- probes$Name
	colnames(M) <- sample_names
	
	U <- U_colour[probes$AddressA, , drop = FALSE]
	rownames(U) <- probes$Name
	colnames(U) <- sample_names
	
	list(M = M, U = U)
}

get_methylation_data <- function(RGset, samples) {
	sample_names <- samples[sampleNames(RGset), "Sample_Name"]
	green <- getGreen(RGset)
	red <- getRed(RGset)

	typeII_signals <- get_signals_type(RGset, "II", green, red, sample_names)
	M <- typeII_signals$M
	U <- typeII_signals$U
	detection_p <- detectionP(RGset, type = "m+u")
	colnames(detection_p) <- as.character(sample_names)

	list(M = M, U = U, detection_p = detection_p)
}

run_workflow <- function(opt) {
	samples <- fread(opt$samples, header=F)

	RGset <- get_RGset(samples, bg_correction = TRUE)
	methylation_data <- get_methylation_data(RGset, samples)
	
	saveRDS(methylation_data, file = opt$output)
}


# ======================= Main

option_list <- list(
        make_option("--samples", type = "character", default = NULL,
                help = "Two-column table of sample names and path to IDAT files.", metavar = "character"),
	make_option(c("-o", "--output"), type = "character", default = "meth_formatted.RDS",
                help = "Output with methylation intensity matrices (RDS format).", metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$samples)) {
        print_help(opt_parser)
        stop("--samples parameter must be provided.", call.=F)
}

run_workflow(opt)

