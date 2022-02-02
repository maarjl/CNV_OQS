# Omics-informed CNV Quality Estimation

## Overview

Recently, studying phenotype associations with copy number variations (CNV) have become more common but the CNV detection step remains challenging and prone to false positive calls. Therefore, additional filtering or quality estimation steps are necessary. This is especially true for CNV detection from widely available genotyping array data. 

Our aim is to develop (and to provide the user with means to develop) an **omics-informed CNV quality score (OQS)** prediction model based on several omics layers from multiple orthogonal sources. Once developed, OQS model can be used to predict CNV quality even in datasets that lack additional omics layers altogether. Its aim is to down-weight CNV calls based on their likelihood of being false positive and, thus, to improve statistical power for follow-up association analyses. 

Our **LINKstudyLINK** focuses on CNV from array data using the popular [PennCNV](http://penncnv.openbioinformatics.org/en/latest/) detection software. We have enclosed our ready-to-use OQS model and the code necessary for the user to directly apply it to their PennCNV data. However, **our quality estimation concepts (and scripts!) are not limited by the choice of detection algorithm**. Even if PennCNV is used, for datasets with at least one additional omics layer available, we strongly encourage the users to build their own custom model tailored for their specific dataset. 

Our workflow (**Figure 1A**) can be summarised by the following steps:

1. **CNV detection** - This step is specific to the choice of your detection method. We used the Hidden Markov Model (HMM)-based [PennCNV](http://penncnv.openbioinformatics.org/en/latest/) software.

2. **CNV quality evaluation based on multiple omics data layers** - Using samples for which additional omics data is available, we calculate an omics-informed quality metric for each CNV. The final metric is combined from up to three individual metrics (with values ranging between 0 and 1) from different omics layers:

	* **whole-genome sequencing (WGS) metric** -- WGS metric can be used to assess the quality of CNVs called from genotyping array intensities if (a subset of) samples also have WGS data available. The metric is calculated as the fraction of CNV region (in basepairs) that can be validated with CNVs called from WGS reads (**Figure 1B**);
	* **gene expression (GE) metric** -- GE metric is based on the assumption that true CNVs alter (deletions decrease and duplications increase) the expression level of genes they overlap, while false calls have no effect on gene expression. The metric captures how extreme a gene expression level of a potential CNV carrier is compared to that of the bulk of the samples assumed to be copy neutral (**Figure 1C**);
	* **methylation intensity (MET) metric** -- MET metric is calculated analogously to GE metric and is based on the assumption that true CNVs alter the overall methylation intensity (sum of methylated and unmethylated intensities) of CpG sites they overlap, while false calls do not (**Figure 1C**).

	If more than one omics layer is used, the final quality metric per CNV is taken to be equal to the 'most extreme' omics-based metric of that CNV (i.e., metric that maximises the `abs(0.5 - metric_value)`).
	
	In case omics data is available for the full set of analysis samples, the combined metrics can be directly applied to the association study. Otherwise, if only a subset of samples are used in omics-based metric calculations, the results should be carried over to the next modelling step. 
	
3. **CNV quality modelling** - This step requires the combined metrics and a set of CNV/sample specific parameters (e.g., CNV length, number of CNVs per sample, etc) output by the CNV detection software, as input. It applies the stepwise model selection approach and chooses the subset of given parameters that best predict the combined CNV quality metric. As a results it builds an omics-informed CNV quality score (OQS) model that can be applied to any set of CNVs and that no longer requires the availability of omics data.

4. **CNV quality prediction using OQS model** - As a final step, OQS model will be used to predict the CNV quality for all analysis samples. The output is a score ranging from 0 to 1, which can be applied to association studies similarly to single nucleotide variant dosage values.

---

![Fig1](figures/fig1_overview.png)

***Figure 1.*** *jou*

---

### How well do the omics-informed quality estimations work?

First, we saw high concordance between our three omics-based metrics. Both GE and MET metrics had high Pearson correlations (R>0.7) when compared to WGS metric. This indicates that either gene expression or methylation intensity data is suitable for CNV quality evaluations if WGS data is not available. All three metrics showed bimodal distribution with modes near 0 and 1 in multiple independent datasets included in our study. Indeed, in all of them >85% of CNVs were evaluated with metric values <0.1 or >0.9. This means that the metrics clearly differentiate between true and false calls for the majority of CNVs.

Secondly, we evaluated the OQS prediction model and demonstrated the improvements achieved using two different approaches: 

* in close family members from three independent datasets, the ‘familial’ CNVs shared between relatives (likely true positives) scored significantly higher with OQS compared to CNVs that were not shared between relatives (mix of true and false positives);
* considering 21 previously published CNV-trait associations (**LINK**) and 89,516 Estonian Biobank samples, we showed that the relative increase in variance explained was up to 34% and 55% when comparing OQS to raw PennCNV output and a previously published CNV quality score (**LINK**), respectively.


Both approaches (with relevant R scripts) are further discussed **LINKbelowLINK**.
For additional information, please see our preprint: **LINK**.

## Workflow setup

The analysis is based on a set of R scripts. All scripts are run from the command line.
R can be downloaded [here](https://www.r-project.org).

The following list of R packages is required for estimating and modelling CNV quality:

~~~r
install.packages("data.table")
~~~

Additionally, these packages are required to run all the code in this repository (including validation analyses presented in our paper):

~~~r
install.packages("something")
~~~

## General usage

### Input

### Calculations of quality metrics based on omics layers

#### Based on WGS

~~~sh
mkdir jou
~~~

#### Based on methylation or gene expression

### Modelling CNV quality

### Scoring CNVs

### PennCNV output conversion

### Applying PennCNV model to your CNV data

## Analyses performed

### Familial CNV analysis

Figure relatives of Omni

![](figures/fig_relatives.png)

### Association analysis per CNV region

Figure of 16p region

![](figures/fig_16p.png)

## Contact and citation

 

