#!/bin/bash
 
#SBATCH -J 1_450_new_CNVD
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -t 192:00:00
#SBATCH --mem=50000

module load java-1.8.0_40
module load tabix-0.2.6
module load samtools-1.2
module load R-3.2.0

export SV_DIR=`cd /path/to/Soft/svtoolkit && pwd`
#SV_TMPDIR=./tmpdir

projectDir=/path/to/full/GenomeStrip
bamList=/path/to/full/GenomeStrip/input/input_batch1_bam_files_rocket.list
GenderMap_file=/path/to/full/GenomeStrip/input/GenderMap_batch1.txt
jobProject=cnvdiscovery_batch1_new
jobQueue=week_batch1_new
mdDir=$projectDir/metadata_batch1
refBundle=/path/to/references/GenomeStrip/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta

#export SV_CLASSPATH="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar"
export SV_CLASSPATH="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
#export LD_LIBRARY_PATH="${SV_DIR}/bwa:${LD_LIBRARY_PATH}"

#logDir=./logs/CNVDiscovery_1_450
#mkdir -p ${logDir} || exit 1

java -cp ${SV_CLASSPATH} \
     org.broadinstitute.gatk.queue.QCommandLine \
     -S ${SV_DIR}/qscript/discovery/cnv/CNVDiscoveryPipeline.q \
     -S ${SV_DIR}/qscript/SVQScript.q \
     -cp ${SV_CLASSPATH} \
     -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
     -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
     -R ${refBundle} \
     -md ${mdDir} \
     -genderMapFile ${GenderMap_file} \
     -runDirectory CNVDiscovery_batch1 \
     -jobProject ${jobProject} \
     -jobQueue ${jobQueue} \
     -jobLogDir CNVDiscovery_batch1/logs \
     -tilingWindowSize 1000 \
     -tilingWindowOverlap 500 \
     -maximumReferenceGapLength 1000 \
     -boundaryPrecision 100 \
     -minimumRefinedLength 500 \
     -I ${bamList} \
   -jobRunner Drmaa \
   -gatkJobRunner Drmaa \
   -jobNative "--mem=8192" \
   -jobNative "--time=192:00" \
   -retry 1 \
     -run \
|| exit 1


