// ATAC-seq pipeline is a sequential pipeline that accomplishes to generate peaks from the fastq data provided

// Here are the steps involved in the ATAC-seq Analysis
// 1. Pre Alignment processing - Trimming the adapters using trim galore.
// 2. Genome alignment using bowtie2.
// 3. Indexing and sorting using bamtools.
// 4. Post Alignment processing - Filtering uninformatice reads: Picard to MarkDuplicates and CollectInsertSizeMetrics.
// 5. Post Alignment processing - Peak calling using macs2.

println " ATAC-seq pipeline"
println "  "
println "Here are the steps involved in the ATAC-seq Analysis"
println "1. Pre Alignment processing - Trimming the adapters using trim galore."
println "2. Genome alignment using bowtie2."
println "3. Indexing and sorting using bamtools."
println "4. Post Alignment processing - Filtering uninformatice reads: Picard to MarkDuplicates and CollectInsertSizeMetrics."
println "5. Post Alignment processing - Peak calling using macs2."


process generate_fastqc_multiqc_reports {
    debug true
    
    input:
    val config_directory
    val fastq_files
    val fastqc_cores

    script:
    """
    echo '${config_directory}/generate_fastqc_reports.sh $fastq_files $fastqc_cores'
    bash ${config_directory}/generate_fastqc_reports.sh $fastq_files $fastqc_cores
    """
}


process trim_galore_adapter_trimming {
    debug true
    
    input:
    val config_directory
    val fastq_files
    val samples_file
    val trim_galore_cores

    script:
    """
    echo '${config_directory}/trim_galore_script.sh $fastq_files $samples_file $trim_galore_cores'
    bash ${config_directory}/trim_galore_script.sh $fastq_files $samples_file $trim_galore_cores
    """
}


process mapping_bowtie2 {
    debug true
    
    input:
    val config_directory
    val fastq_files
    val samples_file

    script:
    """
    echo '${config_directory}/bowtie2_fq_2bam.sh $fastq_files $samples_file'
    bash ${config_directory}/bowtie2_fq_2bam.sh $fastq_files $samples_file
    """
}


process mark_duplicates_picard {
    debug true
    
    input:
    val config_directory
    val fastq_files
    val samples_file
    val picard_filepath

    script:
    """
    echo '${config_directory}/mark_duplicates_bam.sh $fastq_files $samples_file $picard_filepath'
    bash ${config_directory}/mark_duplicates_bam.sh $fastq_files $samples_file $picard_filepath
    """
}


process collect_insert_sizes_picard {
    debug true
    
    input:
    val config_directory
    val fastq_files
    val samples_file
    val picard_filepath

    script:
    """
    echo '${config_directory}/collect_insert_size_metrics.sh $fastq_files $samples_file $picard_filepath'
    bash ${config_directory}/collect_insert_size_metrics.sh $fastq_files $samples_file $picard_filepath
    """
}


process macs2_peak_calling {
    debug true
    
    input:
    val config_directory
    val fastq_files
    val samples_file

    script:
    """
    echo '${config_directory}/macs2_peak_calling.sh $fastq_files $samples_file'
    bash ${config_directory}/macs2_peak_calling.sh $fastq_files $samples_file
    """
}


process generate_bigwig_files {
    debug true
    
    input:
    val config_directory
    val fastq_files
    val hg38gn_filepath

    script:
    """
    echo '${config_directory}/bedGraph_to_bigWig.sh $fastq_files $hg38gn_filepath'
    bash ${config_directory}/bedGraph_to_bigWig.sh $fastq_files $hg38gn_filepath
    """
}




workflow {

    config_directory = params.config_directory
    fastq_files = params.fastq_files
    samples_file = params.samples_file

    // logging config variables
    println " "
    println "Config directory: ${config_directory}"
    println "Fastq files directory: ${fastq_files}"
    println " "

    // Running FastQC and MultiQC Reports
    fastqc_cores = params.fastqc_cores
    println " "
    println "FastQC and MultiQC output directory: ${config_directory}/fastqc_and_multiqc_reports"
    println "No. of cores to be used for generating FastQC reports: ${fastqc_cores}"
    println " "
    generate_fastqc_multiqc_reports(config_directory, fastq_files, fastqc_cores)

    // Running trim-galore
    trim_galore_cores = params.tm_galore_cores
    println " "
    println "Trim Galore output directory: ${config_directory}/trim_galore_output"
    println "No. of cores to be used to trim the adapters: ${trim_galore_cores}"
    println " "
    trim_galore_adapter_trimming(config_directory, fastq_files, samples_file, trim_galore_cores)
    
    // - - - -- - - - - - - ---
    // Note: rename the files -
    // - - - -- - - - - - - ---

    // Running Bowtie2 Mapping
    println " "
    println "Bowtie2 mapping output directory: ${config_directory}/bowtie_output"
    println " "
    mapping_bowtie2(config_directory, fastq_files, samples_file)

    // Running MarkDuplicates (Picard)
    picard_filepath = params.picard_filepath
    println " "
    println "Mark Duplicates using picard output directory: ${config_directory}/mark_duplicate_output"
    println " "
    mark_duplicates_picard(config_directory, fastq_files, samples_file, picard_filepath)
 
    // Running CollectInsertSize (Picard)
    // println " "
    // println "Collect insert size metrics output directory: ${config_directory}/collect_insert_metrics_output"
    // println " "
    // collect_insert_sizes_picard(config_directory, fastq_files, samples_file, picard_filepath)

    // Running MACS2 Peak Calling
    println " "
    println "MACS2 peak-calling output directory: ${config_directory}/macs2_peak_calling_output"
    println " "
    macs2_peak_calling(config_directory, fastq_files, samples_file)

    // Converting .bed to .bedGraph to .bigWig files
    // hg38gn_filepath = params.hg38gn_filepath
    // println " "
    // println "Generating .bigWig files peak-calling output directory: ${config_directory}/macs2_peak_calling_bedgraph_output"
    // println " "
    // generate_bigwig_files(config_directory, fastq_files, hg38gn_filepath)


}
