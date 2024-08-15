#!/bin/bash

set -e

# checking if 4 arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <BAM_INPUT_PATH> <JOBLIST_FQ_FPATH_FNAME> <PICARD_JAR_FPATH> <OUTPUT_BASE_PATH>"
    exit 1
fi

# INPUT PATHS
BAM_INPUT_PATH=$1
JOBLIST_FQ_FPATH_FNAME=$2
PICARD_JAR_FPATH=$3

# OUTPUT PATHS
BOWTIE_INPUT_PATH="$4/bowtie_output"
MARK_DUPLICATES_OUTPUT_PATH="$4/mark_duplicate_output"

# checking if bowtie2-mapping output directory exist
if [ ! -d "$BOWTIE_INPUT_PATH" ]; then
    echo "Error: Directory $BOWTIE_INPUT_PATH does not exist."
    exit 1
fi

mkdir -p $MARK_DUPLICATES_OUTPUT_PATH

# LOG_PATHS
MARK_DUPLICATES_LOG_PATH=${MARK_DUPLICATES_OUTPUT_PATH}/log

mkdir -p $MARK_DUPLICATES_LOG_PATH

for LIBR_NAME in `cat ${JOBLIST_FQ_FPATH_FNAME}`
do
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Running Mark Duplicates ...
    COMMAND="java -jar $PICARD_JAR_FPATH MarkDuplicates INPUT=${BOWTIE_INPUT_PATH}/${LIBR_NAME}_filtered_output.bam OUTPUT=${MARK_DUPLICATES_OUTPUT_PATH}/${LIBR_NAME}_marked_duplicates.bam METRICS_FILE=${MARK_DUPLICATES_OUTPUT_PATH}/${LIBR_NAME}_duplication_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TAGGING_POLICY=All QUIET=true" >> "${MARK_DUPLICATES_LOG_PATH}/${LIBR_NAME}_mapping_step.log"
    echo $COMMAND
    eval $COMMAND
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Mark Duplicates completed.
    echo " "
    echo " "
done
