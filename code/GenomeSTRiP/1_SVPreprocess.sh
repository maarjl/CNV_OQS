#!/bin/bash
 
#SBATCH -J GS
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -t 168:00:00
#SBATCH --mem=50000

module load java-1.8.0_40
module load tabix-0.2.6
module load samtools-1.2
module load R-3.2.0

export SV_DIR=`cd /path/to/Soft/svtoolkit && pwd`
SV_TMPDIR=./tmpdir

projectDir=/path/to/full/GenomeStrip
bamList=/path/to/full/GenomeStrip/input/input_batch1_bam_files_rocket.list
jobProject=preprocess1
jobQueue=week
mdDir=$projectDir/metadata_batch1
refBundle=/path/to/references/GenomeStrip/Homo_sapiens_assembly19/Homo_sapiens_assembly19.fasta

export SV_CLASSPATH="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar"
export LD_LIBRARY_PATH="${SV_DIR}/bwa:${LD_LIBRARY_PATH}"

logDir=./logs/preprocess_batch1
mkdir -p ${logDir} || exit 1


java -cp ${SV_CLASSPATH} \
	org.broadinstitute.gatk.queue.QCommandLine \
	-S ${SV_DIR}/qscript/SVPreprocess.q \
	-S ${SV_DIR}/qscript/SVQScript.q \
	-configFile ${SV_DIR}/conf/genstrip_parameters.txt \
	-tempDir ${SV_TMPDIR} \
	-gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
	-cp ${SV_CLASSPATH} \
	-jobProject ${jobProject} \
	-jobQueue ${jobQueue} \
	-jobLogDir ${logDir} \
	-R ${refBundle} \
	-md ${mdDir} \
	-bamFilesAreDisjoint true \
	-I ${bamList} \
	-jobRunner Drmaa \
	-gatkJobRunner Drmaa \
	-jobNative "--mem=8192" \
	-jobNative "--time=96:00" \
	-P metadata.version:2 \
	-retry 1 \
	-run \
|| exit 1


