#!/bin/bash

set -e

# INPUT PATHS
FQ_INPUT_PATH=$1
JOBLIST_FQ_FPATH_FNAME=$2

# OUTPUT PATHS
BOWTIE_OUTPUT_PATH="$3/bowtie_output"
TRIM_GALORE_OUTPUT_FPATH="$3/trim_galore_output"

# checking if trim-galore output directory exist
if [ ! -d "$TRIM_GALORE_OUTPUT_FPATH" ]; then
    echo "Error: Directory $TRIM_GALORE_OUTPUT_FPATH does not exist."
    exit 1
fi

mkdir -p $BOWTIE_OUTPUT_PATH

NUM_THREADS=45

# INDEX NAMES
IDX_PREFIX=hg38.p13

# FILE NAMES
READ1_SUFFIX="_R1_001"
READ2_SUFFIX="_R2_001"

# LOG_PATHS
BAM_LOG_FPATH=${BOWTIE_OUTPUT_PATH}'/log'

mkdir -p $BAM_LOG_FPATH

for LIBR_NAME in `cat ${JOBLIST_FQ_FPATH_FNAME}`
do
    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Running bowtie2 mapping ...
    READ1_FPATH="${TRIM_GALORE_OUTPUT_FPATH}/$(ls "${TRIM_GALORE_OUTPUT_FPATH}" | grep -i "${LIBR_NAME}${READ1_SUFFIX}" | grep -i 'val_1.fq.gz')"
    READ2_FPATH="${TRIM_GALORE_OUTPUT_FPATH}/$(ls "${TRIM_GALORE_OUTPUT_FPATH}" | grep -i "${LIBR_NAME}${READ2_SUFFIX}" | grep -i 'val_2.fq.gz')"


    # checking if fastq/fq files exist
    if [[ ! -f "$READ1_FPATH" ]]; then
        echo "Error: File $READ1_FPATH does not exist."
        continue
    fi
    if [[ ! -f "$READ2_FPATH" ]]; then
        echo "Error: File $READ2_FPATH does not exist."
        continue
    fi

    echo "bowtie2 --very-sensitive -I 0 -X 1000 -x hg38 -1 \"$READ1_FPATH\" -2 \"$READ2_FPATH\" --fr | samtools view -q 30 -f 2 -F 4 -h - | grep -v 'chrM' | samtools view -bS - > \"${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output.bam\" >> \"${BAM_LOG_FPATH}/${LIBR_NAME}_mapping_step.log\""

    bowtie2 --very-sensitive -I 0 -X 1000 -x hg38 -1 "$READ1_FPATH" -2 "$READ2_FPATH" --fr | samtools view -q 30 -f 2 -F 4 -h - | grep -v 'chrM' | samtools view -bS - > "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output.bam" 2>> "${BAM_LOG_FPATH}/${LIBR_NAME}_mapping_step.log"

    if [ $? -ne 0 ]; then
        echo "Error during bowtie2 mapping for ${LIBR_NAME}"
        continue
    fi

    # sorting, filtering and removing chrM's
    bamtools sort -in "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output.bam" -out "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output_sort.bam"
    bamtools filter -in "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output_sort.bam" -out "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output_filt.bam" -mapQuality ">=30" -isProperPair true
    samtools view -h "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output_filt.bam" | grep -v "chrM" | samtools view -b -o "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output.bam"
    rm "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output_sort.bam" "${BOWTIE_OUTPUT_PATH}/${LIBR_NAME}_filtered_output_filt.bam"

    echo `date '+%F %H:%M:%S'` ${LIBR_NAME} Bowtie2 mapping completed.
    echo " "
    echo " "
done
