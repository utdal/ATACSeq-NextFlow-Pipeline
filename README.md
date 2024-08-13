# ATACSeq Data Processing Pipeline using NextFlow
This is an automated workflow pipeline for analyzing and processing ATAC-seq data, implemented primarily in bash, and wrapped in a NextFlow workflow to characterize the chromatin landscape in bulk ATAC-seq samples. Here are the steps for data processing:
1. [Completed] Running Trim galore to cut the adapters
2. [Completed] Running alignment to the reference genome using Bowtie2
3. [Completed] Running filtering using Samtools
4. [Completed] Running mark duplicates using picard
5. [Completed] Running peak calling using MACS2
6. [In-progress] Calculating TSSe score
7. [In-progress] Generating bigWig and heatmap using Deeptools

![ATACSeq NextFlow Pipeline](misc/ATACSeqpipeline.png)

This tool is used to process bulk ATAC-seq data by mapping paired-end reads to a reference genome and identifying areas of open chromatin after peak calling. This tool generates files that can be visualized on a genome browser.

Running the tool is pretty straight forward, however a good understanding of `bash` is recommended. Please familiarize yourself with data types, basic functions, data structures in each language.

## Installation/Setup of Visium NextFlow Pipeline:
You can install Visium NextFlow Pipeline via git:
```
git clone https://github.com/utdal/ATACSeq-NextFlow-Pipeline
```

To execute the tool, essential modifications need to be made to the file(s):
```
a) pipeline.config
b) atac_seq_samples.txt
```

> Note:
> 1. Install [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt), [MultiQC](https://multiqc.info/docs/getting_started/installation/), [cutadapt](https://cutadapt.readthedocs.io/en/v3.7/installation.html#installation-with-conda), [trim-galore](https://github.com/FelixKrueger/TrimGalore), and macs2 packages.
> 2. To run `MarkDuplicates`, you will need the Picard Java archive file `picard.jar`, which can be downloaded from the [Broad Institute's](https://github.com/broadinstitute/picard/releases/tag/3.2.0) website. Make sure to update the file path to this archive in the `pipeline.config` file.
> 3. Before executing the pipeline, you must build the **Bowtie2 index** from the reference genome and place it in the config directory: `params.config_directory = '/path/to/config'`.
> Download the reference genome: `hg38canon.fa` and, to build the index, execute: `bowtie2-build hg38canon.fa hg38`
> 
> Here **hg38** is Bowtie2 human genome index directory, that needs to be updated in the `pipeline.config` file.

#### Running the Tool:
Here is an example of how to run the pipeline:
1. Command to run the pipeline:
   ```
   nextflow run atac_seq_nextflow_pipeline.nf -c pipeline.config
   ```
2. Command to re-run from fail-point:
   ```
   nextflow run atac_seq_nextflow_pipeline.nf -c pipeline.config -resume
   ```

The results generated are stored in the `params.config_directory = '/path/to/config'` directory, as mentioned in the `pipeline.config` file.
