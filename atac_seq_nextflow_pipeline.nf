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
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val fastq_files
    val fastqc_cores

    output:
    val "${config_directory}/fastqc_and_multiqc_reports", emit: fastqc_output

    script:
    """
    echo '${config_directory}/generate_fastqc_reports.sh $fastq_files $fastqc_cores $config_directory'
    bash ${config_directory}/generate_fastqc_reports.sh $fastq_files $fastqc_cores $config_directory
    """
}


process trim_galore_adapter_trimming {
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val fastq_files
    val samples_file
    val trim_galore_cores
    val fastqc_output

    output:
    val "${config_directory}/trim_galore_output", emit: trimmed_files

    script:
    """
    echo '${config_directory}/trim_galore_script.sh $fastq_files $samples_file $trim_galore_cores $config_directory'
    bash ${config_directory}/trim_galore_script.sh $fastq_files $samples_file $trim_galore_cores $config_directory
    """
}


process mapping_bowtie2 {
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val samples_file
    val genome_index_directory
    val trimmed_files

    output:
    val "${config_directory}/bowtie_output", emit: mapped_files

    script:
    """
    echo '${config_directory}/bowtie2_fq_2bam.sh $trimmed_files $samples_file $config_directory $genome_index_directory'
    bash ${config_directory}/bowtie2_fq_2bam.sh $trimmed_files $samples_file $config_directory $genome_index_directory
    """
}


process mark_duplicates_picard {
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val samples_file
    val picard_filepath
    val mapped_files

    output:
    val "${config_directory}/mark_duplicate_output", emit: marked_dup_files

    script:
    """
    echo '${config_directory}/mark_duplicates_bam.sh $mapped_files $samples_file $picard_filepath $config_directory'
    bash ${config_directory}/mark_duplicates_bam.sh $mapped_files $samples_file $picard_filepath $config_directory
    """
}


process collect_insert_sizes_picard {
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val samples_file
    val picard_filepath
    val marked_dup_files

    output:
    val "${config_directory}/collect_insert_metrics_output", emit: insert_size_metrics

    script:
    """
    echo '${config_directory}/collect_insert_size_metrics.sh $marked_dup_files $samples_file $picard_filepath $config_directory'
    bash ${config_directory}/collect_insert_size_metrics.sh $marked_dup_files $samples_file $picard_filepath $config_directory
    """
}


process macs2_peak_calling {
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val samples_file
    val marked_dup_files

    output:
    val "${config_directory}/macs2_peak_calling_output", emit: peak_files

    script:
    """
    echo '${config_directory}/macs2_peak_calling.sh $marked_dup_files $samples_file $config_directory'
    bash ${config_directory}/macs2_peak_calling.sh $marked_dup_files $samples_file $config_directory
    """
}


process generate_bigwig_files {
    debug true
    errorStrategy 'terminate'
    
    input:
    val config_directory
    val hg38gn_filepath
    val peak_files

    output:
    val "${config_directory}/macs2_peak_calling_bedgraph_output", emit: bigwig_files

    script:
    """
    echo '${config_directory}/bedGraph_to_bigWig.sh $peak_files $hg38gn_filepath $config_directory'
    bash ${config_directory}/bedGraph_to_bigWig.sh $peak_files $hg38gn_filepath $config_directory
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
    fastqc_output = generate_fastqc_multiqc_reports(config_directory, fastq_files, fastqc_cores)

    // Running trim-galore
    trim_galore_cores = params.tm_galore_cores
    println " "
    println "Trim Galore output directory: ${config_directory}/trim_galore_output"
    println "No. of cores to be used to trim the adapters: ${trim_galore_cores}"
    println " "
    trimmed_files = trim_galore_adapter_trimming(config_directory, fastq_files, samples_file, trim_galore_cores, fastqc_output.fastqc_output)
    
    // Running Bowtie2 Mapping
    genome_index_directory = params.genome_index_directory
    println " "
    println "Bowtie2 mapping output directory: ${config_directory}/bowtie_output"
    println " "
    mapped_files = mapping_bowtie2(config_directory, samples_file, genome_index_directory, trimmed_files.trimmed_files)

    // Running MarkDuplicates (Picard)
    picard_filepath = params.picard_filepath
    println " "
    println "Mark Duplicates using picard output directory: ${config_directory}/mark_duplicate_output"
    println " "
    marked_dup_files = mark_duplicates_picard(config_directory, samples_file, picard_filepath, mapped_files.mapped_files)
 
    // Running CollectInsertSize (Picard)
    println " "
    println "Collect insert size metrics output directory: ${config_directory}/collect_insert_metrics_output"
    println " "
    insert_size_metrics = collect_insert_sizes_picard(config_directory, samples_file, picard_filepath, marked_dup_files.marked_dup_files)

    // Running MACS2 Peak Calling
    println " "
    println "MACS2 peak-calling output directory: ${config_directory}/macs2_peak_calling_output"
    println " "
    peak_files = macs2_peak_calling(config_directory, samples_file, marked_dup_files.marked_dup_files)

    // Converting .bed to .bedGraph to .bigWig files
    hg38gn_filepath = params.hg38gn_filepath
    println " "
    println "Generating .bigWig files peak-calling output directory: ${config_directory}/macs2_peak_calling_bedgraph_output"
    println " "
    bigwig_files = generate_bigwig_files(config_directory, hg38gn_filepath, peak_files.peak_files)

}
