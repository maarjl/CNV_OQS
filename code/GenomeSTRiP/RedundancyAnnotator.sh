#!/bin/bash
 
#SBATCH -J RA
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -t 10:00:00
#SBATCH --mem=50000

module load java-1.8.0_40
#module load tabix-0.2.6
#module load samtools-1.2
module load R-3.2.0

export SV_DIR=`cd /path/to/Soft/svtoolkit && pwd`
export SV_CLASSPATH="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"


projectDir=/path/to/full/GenomeStrip
mdDir=$projectDir/metadata_1_2284
refBundle=/path/to/references/GenomeStrip/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta

input_vcf=$projectDir/CNVDiscovery/ALL/merge/RundancyAnnotator/genotypes_1_2284.vcf
output_vcf=$projectDir/CNVDiscovery/ALL/merge/RundancyAnnotator/genotypes_1_2284_RA.vcf
reportDir=$projectDir/CNVDiscovery/ALL/merge/RundancyAnnotator

java -cp ${SV_CLASSPATH} \
     org.broadinstitute.sv.main.SVAnnotator \
     -A Redundancy \
     -R ${refBundle} \
     -vcf ${input_vcf} \
     -comparisonFile ${input_vcf} \
     -duplicateOverlapThreshold 0.5 \
     -O ${output_vcf} \
     -writeReport true \
     -reportDirectory ${reportDir}


