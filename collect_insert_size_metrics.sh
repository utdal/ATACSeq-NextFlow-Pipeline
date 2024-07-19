#!/bin/bash

# INPUT PATHS
BAM_INPUT_PATH=$1
JOBLIST_FQ_FPATH_FNAME=$2
PICARD_JAR_FPATH=$3

# OUTPUT PATHS
MARK_DUPLICATES_OUTPUT_PATH="$BAM_INPUT_PATH/mark_duplicate_output"
COLLECT_INSERT_SIZES_OUTPUT_PATH="$BAM_INPUT_PATH/collect_insert_metrics_output"

mkdir -p $COLLECT_INSERT_SIZES_OUTPUT_PATH


# LOG_PATHS
COLLECT_INSERT_SIZES_LOG_PATH=${COLLECT_INSERT_SIZES_OUTPUT_PATH}'/log'

mkdir -p $COLLECT_INSERT_SIZES_LOG_PATH

for LIBR_NAME in `cat ${JOBLIST_FQ_FPATH_FNAME}`
do
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Running Mark Duplicates ...
    COMMAND="java -jar $PICARD_JAR_FPATH CollectInsertSizeMetrics INPUT=${MARK_DUPLICATES_OUTPUT_PATH}/${LIBR_NAME}_marked_duplicates.bam OUTPUT=${COLLECT_INSERT_SIZES_OUTPUT_PATH}/${LIBR_NAME}_insert_size_metric.tsv H=${COLLECT_INSERT_SIZES_OUTPUT_PATH}/${LIBR_NAME}_insert_size_metric.pdf" >> "${COLLECT_INSERT_SIZES_LOG_PATH}/${LIBR_NAME}_collect_insert_metrics_step.log"
    echo $COMMAND
    eval $COMMAND
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Mark Duplicates completed.
    echo " "
    echo " "
done