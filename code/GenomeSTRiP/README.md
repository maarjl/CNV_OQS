[Genome STRiP](https://software.broadinstitute.org/software/genomestrip/) is a pipeline for structural variant (SV), including CNVs, detection from sequencing reads. Here we present an overview of the workflow used on 2,284 Estonian whole-genome sequencing (WGS) samples. First steps are executed in five batches of maximum 450 samples each. The complete tutorial can be found [here](https://software.broadinstitute.org/software/genomestrip/sites/default/files/materials/GATKWorkshop_GenomeSTRiP_tutorial_July2013.pdf). Besides Genome STRiP, also GATK and bcftools are required to run the full pipeline.

1. `1_SVPreprocess.sh`: Preprocessing of BAM files (per batch)
2. `2_SVDiscovery.sh`: Discovery of SV regions (per batch)
3. `3_SVGenotypes.sh`: Genotyping of regions found in previous step (per batch)
4. Merging of batch-specific `3_SVGenotypes.sh` output files using `bcftools`
5. `RedundancyAnnotator.sh`: Removal of duplicates based on site overlap greater than 50%
6. Removal of sites with duplicate score >0: `
bcftools view --exclude 'GSDUPLICATESCORE>0' genotypes_1_2284_RA.vcf -Ov -o EGCUT_CNVs_final.vcf
`
7. `VariantFiltration.sh`: Filter variants based on various cutoff values 
