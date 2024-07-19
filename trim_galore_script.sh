#!/bin/bash


# INPUT PATHS
FASTQ_INPUT_PATH=$1
JOBLIST_FPATH_FNAME=$2

# OUTPUT PATHS
TRIM_GALORE_OUTPUT_FPATH="$1/trim_galore_output"
MAP_OUTPUT_FPATH="$1/bam_output"

mkdir -p $TRIM_GALORE_OUTPUT_FPATH $MAP_OUTPUT_FPATH

NUM_THREADS=$3

# INDEX NAMES
IDX_PREFIX=hg38.p13

# FILE NAMES
READ1_SUFFIX="_R1_001"
READ2_SUFFIX="_R2_001"

# LOG_PATHS
TRIM_GALORE_LOG_FPATH=${TRIM_GALORE_OUTPUT_FPATH}'/log'
BAM_LOG_FPATH=${MAP_OUTPUT_FPATH}'/log'
PEAKS_LOG_FPATH=${PKCALL_OUTPUT_FPATH}/'log'

mkdir -p $TRIM_GALORE_LOG_FPATH $BAM_LOG_FPATH

for LIBR_NAME in `cat ${JOBLIST_FPATH_FNAME}`
do
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Running trim-galore ...

    #( trim_galore --paired --cores $NUM_THREADS --output $TRIM_GALORE_OUTPUT_FPATH ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ1_SUFFIX}.fastq.gz ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ2_SUFFIX}.fastq.gz ) 2> ${TRIM_GALORE_LOG_FPATH}/${LIBR_NAME}_mapping_step.log
    #echo '"trim_galore --paired --cores $NUM_THREADS --output $TRIM_GALORE_OUTPUT_FPATH ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ1_SUFFIX}.fastq.gz ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ2_SUFFIX}.fastq.gz" >> "${TRIM_GALORE_LOG_FPATH}/${LIBR_NAME}_mapping_step.log"'
    
    ( trim_galore --paired --cores $NUM_THREADS --output $TRIM_GALORE_OUTPUT_FPATH ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ1_SUFFIX}.fastq.gz ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ2_SUFFIX}.fastq.gz ) 2> ${TRIM_GALORE_LOG_FPATH}/${LIBR_NAME}_mapping_step.log

    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Trim-galore completed.
done