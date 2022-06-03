#!/bin/bash
 
#SBATCH -J SVG_1_450
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -t 168:00:00
#SBATCH --mem=5000

module load java-1.8.0_40
module load R-3.2.0

export SV_DIR=`cd /path/to/Soft/svtoolkit && pwd`
export SV_CLASSPATH="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"


projectDir=/path/to/full/GenomeStrip
mdDir=$projectDir/metadata_batch1
refBundle=/path/to/references/GenomeStrip/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta
bamList=$projectDir/CNVDiscovery/ALL/merge/pass_samples/pass_batch1_bam_files_rocket.list
GenderMap_file=$projectDir/input/GenderMap_batch1.txt

input_vcf=$projectDir/CNVDiscovery/ALL/merge/unique_genotypes_col18_1_2284.vcf

java -Xmx4g -cp ${SV_CLASSPATH} \
     org.broadinstitute.gatk.queue.QCommandLine \
     -S ${SV_DIR}/qscript/SVGenotyper.q \
     -S ${SV_DIR}/qscript/SVQScript.q \
     -cp ${SV_CLASSPATH} \
     -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
     -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
     -R ${refBundle} \
     -I ${bamList} \
     -genderMapFile ${GenderMap_file} \
     -md ${mdDir} \
     -runDirectory SVGenotyper_batch1 \
     -jobLogDir SVGenotyper_batch1/logs \
     -vcf ${input_vcf} \
     -O run1/genotypes_batch1.vcf \
     -parallelJobs 450 \
   -jobRunner Drmaa \
   -gatkJobRunner Drmaa \
   -jobNative "--mem=3000" \
   -jobNative "--time=168:00" \
     -run
