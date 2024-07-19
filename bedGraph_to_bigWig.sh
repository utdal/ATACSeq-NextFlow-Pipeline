#!/bin/bash

# INPUT PATHS
FASTQ_INPUT_PATH=$1
HG38_CHROM_SIZES_PATH=$2

# OUTPUT PATHS
MACS2_OUTPUT_PATH="$FASTQ_INPUT_PATH/macs2_peak_calling_output"
MACS2_BEDGRAPH_OUTPUT_PATH="$FASTQ_INPUT_PATH/macs2_peak_calling_bedgraph_output"
BIGWIG_OUTPUT_PATH="$FASTQ_INPUT_PATH/bigwig_output"

mkdir -p $MACS2_BEDGRAPH_OUTPUT_PATH $BIGWIG_OUTPUT_PATH

find "$MACS2_OUTPUT_PATH" -type f -name "*.bed" -exec cp {} "$MACS2_BEDGRAPH_OUTPUT_PATH" \;

# LOG_PATHS
BIGWIG_LOG_PATH=${BIGWIG_OUTPUT_PATH}/log
mkdir -p $BIGWIG_LOG_PATH

for LIBR_NAME in "$MACS2_BEDGRAPH_OUTPUT_PATH"/*.bed
do
    BASE_FILE_NAME=$(basename "$LIBR_NAME" .bed)
    echo " "
    echo `date '+%F %H:%M:%S'` Converting $BASE_FILE_NAME.bed to $BASE_FILE_NAME.bigWig Calling ...
    bedtools genomecov -bg -i $LIBR_NAME -g $HG38_CHROM_SIZES_PATH > $MACS2_BEDGRAPH_OUTPUT_PATH/$BASE_FILE_NAME.bedgraph
    COMMAND="bedGraphToBigWig $MACS2_BEDGRAPH_OUTPUT_PATH/$BASE_FILE_NAME.bedgraph $HG38_CHROM_SIZES_PATH $BIGWIG_OUTPUT_PATH/$BASE_FILE_NAME.bigWig" >> "${BIGWIG_LOG_PATH}/${BASE_FILE_NAME}_conv_bedgraph_to_bigwig_step.log"
    echo $COMMAND
    eval $COMMAND
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Converting $LIBR_NAME to $BASE_FILE_NAME.bigWig completed.
    echo " "
    echo " "
done