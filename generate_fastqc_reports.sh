#!/bin/bash

set -e

NUM_THREADS=$2 
FASTQ_INPUT_PATH=$1
output_dir=$FASTQ_INPUT_PATH/"fastqc_and_multiqc_reports"

if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

cd $FASTQ_INPUT_PATH


for file in $FASTQ_INPUT_PATH/*.fastq.gz; do
    echo "fastqc -t $NUM_THREADS $file --outdir='$output_dir'"
    fastqc -t $NUM_THREADS $file --outdir=$output_dir
done

echo "Successfully completed running fastqc reports on $FASTQ_INPUT_PATH samples."
echo " "
echo " "


cd $output_dir

multiqc .

echo "Successfully completed running multiqc reports on $FASTQ_INPUT_PATH samples."
echo " "
echo " "

echo "FastQC and MultiQC reports are stored here: $output_dir."
