#!/bin/bash

# INPUT PATHS
BAM_INPUT_PATH=$1
JOBLIST_FQ_FPATH_FNAME=$2

# OUTPUT PATHS
MARK_DUPLICATES_OUTPUT_PATH="$BAM_INPUT_PATH/mark_duplicate_output"
MACS2_OUTPUT_PATH="$BAM_INPUT_PATH/macs2_peak_calling_output"

mkdir -p $MACS2_OUTPUT_PATH

# LOG_PATHS
MACS2_LOG_PATH=${MACS2_OUTPUT_PATH}/log

mkdir -p $MACS2_LOG_PATH

for LIBR_NAME in `cat ${JOBLIST_FQ_FPATH_FNAME}`
do
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Running MACS2 Peak Calling ...
    mkdir -p ${MACS2_OUTPUT_PATH}/${LIBR_NAME}
    cd ${MACS2_OUTPUT_PATH}/${LIBR_NAME}
    COMMAND="macs2 callpeak -t ${MARK_DUPLICATES_OUTPUT_PATH}/${LIBR_NAME}_marked_duplicates.bam -n $LIBR_NAME -f BAMPE --nomodel --shift -100 --extsize 200 --keep-dup all -g hs --call-summits --d-min 20 --qvalue '0.05' --bdg" >> "${MACS2_LOG_PATH}/${LIBR_NAME}_mapping_step.log"
    echo $COMMAND
    eval $COMMAND
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} MACS2 Peak Calling completed.
    echo " "
    echo " "
done