#!/bin/bash

set -e

# INPUT PATHS
FASTQ_INPUT_PATH=$1
JOBLIST_FPATH_FNAME=$2

# OUTPUT PATHS
TRIM_GALORE_OUTPUT_FPATH="$4/trim_galore_output"
mkdir -p $TRIM_GALORE_OUTPUT_FPATH

NUM_THREADS=$3

# FILE NAMES
READ1_SUFFIX="_R1_001"
READ2_SUFFIX="_R2_001"

# LOG_PATHS
TRIM_GALORE_LOG_FPATH=${TRIM_GALORE_OUTPUT_FPATH}'/log'
mkdir -p $TRIM_GALORE_LOG_FPATH

for LIBR_NAME in `cat ${JOBLIST_FPATH_FNAME}`
do
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Running trim-galore ...

    #( trim_galore --paired --cores $NUM_THREADS --output $TRIM_GALORE_OUTPUT_FPATH ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ1_SUFFIX}.fastq.gz ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ2_SUFFIX}.fastq.gz ) 2> ${TRIM_GALORE_LOG_FPATH}/${LIBR_NAME}_mapping_step.log
    echo '"trim_galore --paired --cores $NUM_THREADS --output $TRIM_GALORE_OUTPUT_FPATH ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ1_SUFFIX}.fastq.gz ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ2_SUFFIX}.fastq.gz" >> "${TRIM_GALORE_LOG_FPATH}/${LIBR_NAME}_mapping_step.log"'
    ( trim_galore --paired --cores $NUM_THREADS --output $TRIM_GALORE_OUTPUT_FPATH ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ1_SUFFIX}.fastq.gz ${FASTQ_INPUT_PATH}/${LIBR_NAME}${READ2_SUFFIX}.fastq.gz ) 2> ${TRIM_GALORE_LOG_FPATH}/${LIBR_NAME}_mapping_step.log

    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Trim-galore completed.
done
