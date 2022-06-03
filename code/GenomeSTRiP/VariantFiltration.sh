#!/bin/bash
 
#SBATCH -J plot
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -t 1:00:00
#SBATCH --mem=20000

module load java-1.8.0_40

ref_fasta=/path/to/references/GenomeStrip/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta
dir=/path/to/full/GenomeStrip/CNVDiscovery/ALL/merge/RundancyAnnotator/summary
soft=/path/to/software/nextgen_pipeline/programs


java -Xmx10g -jar $soft/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $ref_fasta \
--variant $dir/EGCUT_CNVs_final.vcf \
--filterExpression "GSELENGTH < 1000" \
--filterName "ALIGN" \
--filterExpression "GSCALLRATE < 0.9" \
--filterName "CALLRATE" \
--filterExpression "GSCLUSTERSEP <= 5.0" \
--filterName "CLUSTERSEP" \
--filterExpression "GSM1 <= 0.5 || GSM1 >= 2.0" \
--filterName "GTDEPTH" \
--filterExpression "(1.0 * GSELENGTH / GCLENGTH) < 0.5" \
--filterName "DENSITY" \
--filterExpression "GSVDJFRACTION > 0.0" \
--filterName "VDJREGION" \
-o $dir/EGCUT_CNVs_filtered_mild.vcf \
-log $dir/log.GATK.VariantFiltration.txt




