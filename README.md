# ATACSeq Data Processing Pipeline using NextFlow
This is an automated workflow pipeline for analyzing and processing ATAC-seq data, implemented primarily in bash, and wrapped in a NextFlow workflow to generate peak-calling and TSSe Score calculation. Here are the steps for data processing:
1. [In-progress] Running Trim galore to cut the adapters
2. [In-progress] Running alignment to the reference genome using Bowtie2
3. [In-progress] Running filtering using Samtools
4. [In-progress] Running mark duplicates using picard
5. [In-progress] Running peak calling using MACS2 and TSSe score calculation
6. [In-progress] Generating heatmap using Deeptools

![ATACSeq NextFlow Pipeline](misc/ATACSeqpipeline.png)

This tool is predominantly used for the analysis of ATAC-seq data to obtain TSSe Scores and various plots.

Running the tool is pretty straight forward, however a good understanding of `bash` is recommended. Please familiarize yourself with data types, basic functions, data structures in each language.

## Installation/Setup of Visium NextFlow Pipeline:
You can install Visium NextFlow Pipeline via git:
```
git clone https://github.com/utdal/ATACSeq-NextFlow-Pipeline
```

To execute the tool, essential modifications need to be made to the file(s):
```
> pipeline.config
```


#### Running the Tool
Here is an example of how to run the pipeline:
1. Command to run the pipeline:
   ```
   nextflow run atacseq_nextflow_pipeline.nf -c pipeline.config
   ```
2. Command to re-run from fail-point:
   ```
   nextflow run atacseq_nextflow_pipeline.nf -c pipeline.config -resume
   ```
