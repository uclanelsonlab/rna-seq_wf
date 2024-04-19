version 1.0

task picard_CollectMultipleMetrics {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Collect multiple classes of metrics. This 'meta-metrics' tool runs one or more of the metrics collection modules at the same time to cut down on the time spent reading in data from input files. https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard-"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam_file
        File bam_index
        String prefix
        String? java_mem
        String? picard_docker
    }
    String actual_java_mem=select_first([java_mem, "-Xmx100g"])
    String actual_picard_docker=select_first([picard_docker, "broadinstitute/picard:latest"])
    command {
        set -uexo pipefail
        java "${actual_java_mem}" -jar /usr/picard/picard.jar CollectMultipleMetrics \
        -I ~{bam_file} \
        -R ~{fasta} \
        -O ~{prefix}_multiple_metrics
    }
    runtime {
        docker: actual_picard_docker
        dx_instance_type: "mem2_ssd1_v2_x32"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "4H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File alignment_summary_metrics = "${prefix}_multiple_metrics.alignment_summary_metrics"
        File base_distribution_by_cycle = "${prefix}_multiple_metrics.base_distribution_by_cycle.pdf"
        File base_distribution_by_cycle_metrics = "${prefix}_multiple_metrics.base_distribution_by_cycle_metrics"
        File insert_size_histogram = "${prefix}_multiple_metrics.insert_size_histogram.pdf"
        File insert_size_metrics = "${prefix}_multiple_metrics.insert_size_metrics"
        File quality_by_cycle = "${prefix}_multiple_metrics.quality_by_cycle.pdf"
        File quality_by_cycle_metrics = "${prefix}_multiple_metrics.quality_by_cycle_metrics"
        File quality_distribution = "${prefix}_multiple_metrics.quality_distribution.pdf"
        File quality_distribution_metrics = "${prefix}_multiple_metrics.quality_distribution_metrics"
        File read_length_histogram = "${prefix}_multiple_metrics.read_length_histogram.pdf"
    }
}

task picard_CollectWgsMetrics {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined. https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard-"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam_file
        File bam_index
        String prefix
        String? java_mem
        String? picard_docker
    }
    String actual_java_mem=select_first([java_mem, "-Xmx100g"])
    String actual_picard_docker=select_first([picard_docker, "broadinstitute/picard:latest"])
    command {
        set -uexo pipefail
        java "${actual_java_mem}" -jar /usr/picard/picard.jar CollectWgsMetrics \
            -I ~{bam_file} \
            -R ~{fasta} \
            -O ~{prefix}_collect_wgs_metrics
        ls *
        ls ~{prefix}_collect_wgs_metrics*
    }
    runtime {
        docker: actual_picard_docker
        dx_instance_type: "mem2_ssd1_v2_x32"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "6H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File collect_wgs_metrics = "${prefix}_collect_wgs_metrics"
    }
}

workflow picard_qc {
    parameter_meta {
        fasta: {
            description: "Unziped reference FASTA file stored by DNAnexus",
            extension: ".fasta"
        }
        fasta_fai: {
            description: "Index for reference FASTA file stored by DNAnexus",
            extension: ".fai"
        }
        fasta_dict: {
            description: "Dictionary for reference FASTA file stored by DNAnexus",
            extension: ".dict"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        bam_file: {
            description: "Last BAM file created by the pipeline",
            extension: ".bam"
        }
        bam_index: {
            description: "Index for the BAM file used",
            extension: ".bai"
        }
        picard_docker: {
            description: "TAR zip docker image from Picard stored on DNAnexus ()",
            extension: ".tar.gz"
        }
        java_mem: {
            description: "String for Java memory (default: -Xmx100g)"
        }
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam_file
        File bam_index
        String prefix
        String? java_mem
        String? picard_docker
    }
    call picard_CollectMultipleMetrics {
        input:
            fasta = fasta,
            fasta_fai = fasta_fai,
            fasta_dict = fasta_dict,
            bam_file = bam_file,
            bam_index = bam_index,
            prefix = prefix,
            java_mem = java_mem,
            picard_docker = picard_docker
    }
    call picard_CollectWgsMetrics {
        input:
            fasta = fasta,
            fasta_fai = fasta_fai,
            fasta_dict = fasta_dict,
            bam_file = bam_file,
            bam_index = bam_index,
            prefix = prefix,
            java_mem = java_mem,
            picard_docker = picard_docker
    }
    output {
        File alignment_summary_metrics = picard_CollectMultipleMetrics.alignment_summary_metrics
        File base_distribution_by_cycle = picard_CollectMultipleMetrics.base_distribution_by_cycle
        File base_distribution_by_cycle_metrics = picard_CollectMultipleMetrics.base_distribution_by_cycle_metrics
        File insert_size_histogram = picard_CollectMultipleMetrics.insert_size_histogram
        File insert_size_metrics = picard_CollectMultipleMetrics.insert_size_metrics
        File quality_by_cycle = picard_CollectMultipleMetrics.quality_by_cycle
        File quality_by_cycle_metrics = picard_CollectMultipleMetrics.quality_by_cycle_metrics
        File quality_distribution = picard_CollectMultipleMetrics.quality_distribution
        File quality_distribution_metrics = picard_CollectMultipleMetrics.quality_distribution_metrics
        File read_length_histogram = picard_CollectMultipleMetrics.read_length_histogram
        File collect_wgs_metrics = picard_CollectWgsMetrics.collect_wgs_metrics
    }
}