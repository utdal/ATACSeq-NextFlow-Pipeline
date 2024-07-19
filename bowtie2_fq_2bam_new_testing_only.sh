#!/bin/bash

# INPUT PATHS
FQ_INPUT_PATH=$1
JOBLIST_FQ_FPATH_FNAME=$2

# OUTPUT PATHS
BOWTIE_OUTPUT_PATH="$FQ_INPUT_PATH/bowtie_output"
TRIM_GALORE_OUTPUT_FPATH="$FQ_INPUT_PATH/trim_galore_output"

echo "${BOWTIE_OUTPUT_PATH} ${TRIM_GALORE_OUTPUT_FPATH}"

mkdir -p $BOWTIE_OUTPUT_PATH $TRIM_GALORE_OUTPUT_FPATH

NUM_THREADS=5

# INDEX NAMES
IDX_PREFIX=hg38.p13

# FILE NAMES
READ1_SUFFIX="_R1_001"
READ2_SUFFIX="_R2_001"

# LOG_PATHS
TRIM_GALORE_LOG_FPATH=${TRIM_GALORE_OUTPUT_FPATH}'/log'
BAM_LOG_FPATH=${BOWTIE_OUTPUT_PATH}'/log'
PEAKS_LOG_FPATH=${PKCALL_OUTPUT_FPATH}/'log'

echo "${TRIM_GALORE_LOG_FPATH} ${BAM_LOG_FPATH}"
mkdir -p $TRIM_GALORE_LOG_FPATH $BAM_LOG_FPATH

for LIBR_NAME in `cat ${JOBLIST_FQ_FPATH_FNAME}`
do
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Running trim-galore ...
    echo "bowtie2 --very-sensitive -I 0 -X 1000 -x hg38 -1 \"${TRIM_GALORE_OUTPUT_FPATH}/${LIBR_NAME}${READ1_SUFFIX}.fq.gz\" -2 \"${TRIM_GALORE_OUTPUT_FPATH}/${LIBR_NAME}${READ2_SUFFIX}.fq.gz\" --fr | samtools view -q 30 -f 2 -F 4 -h - | grep -v 'chrM' | samtools view -bS - > \"${LIBR_NAME}_filtered_output.bam\" >> \"${BAM_LOG_FPATH}/${LIBR_NAME}_mapping_step.log\""
    bowtie2 --very-sensitive -I 0 -X 1000 -x hg38 -1 "${TRIM_GALORE_OUTPUT_FPATH}/${LIBR_NAME}${READ1_SUFFIX}.fq.gz" -2 "${TRIM_GALORE_OUTPUT_FPATH}/${LIBR_NAME}${READ2_SUFFIX}.fq.gz" --fr | samtools view -q 30 -f 2 -F 4 -h - | grep -v 'chrM' | samtools view -bS - > "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output.bam" 2>> "${BAM_LOG_FPATH}/${LIBR_NAME}_mapping_step.log"
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Trim-galore completed.
    echo " "
    echo " "
done

echo "The bowtie output is stored here: "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output.bam""

