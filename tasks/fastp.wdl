version 1.0

task fastp_stats {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Fastp gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (A/C/G/T). You get information about adapter contamination and other overrepresented sequences."
    }
    input {
        File fastq_r1
        File fastq_r2
        String prefix
        Int? threads_n
        String? fastp_docker
    }

    Int actual_threads_n=select_first([threads_n, 12])
    String actual_fastp_docker=select_first([fastp_docker, "quay.io/biocontainers/fastp:0.23.3--h5f740d0_0"])
    command {
        set -uexo pipefail
        fastp --thread ~{actual_threads_n} \
        -i ~{fastq_r1} \
        -I ~{fastq_r2} \
        -o ~{prefix}_fastp_R1.fastq.gz \
        -O ~{prefix}_fastp_R2.fastq.gz \
        --html ~{prefix}_fastp.html \
        --json ~{prefix}_fastp.json
    }
    runtime {
        docker: actual_fastp_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "5H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File trimm_r1 = "${prefix}_fastp_R1.fastq.gz"
        File trimm_r2 = "${prefix}_fastp_R2.fastq.gz"
        File fastp_html = "${prefix}_fastp.html"
        File fastp_json = "${prefix}_fastp.json"
    }
}

workflow fastp_stats_wf {
    parameter_meta {
        fastq_r1: {
            description: "Sample R1 FASTQ file",
            extension: ".R1.fastq.gz"
        }
        fastq_r2: {
            description: "Sample R2 FASTQ file",
            extension: ".R2.fastq.gz"
        }
        fastp_docker: {
            description: "TAR zip docker image from fastp stored on DNAnexus (file-GVzvBQ802k8k4zPYKqpfz1vZ)",
            extension: ".tar.gz"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
    }
    input {
        File fastq_r1
        File fastq_r2
        String prefix
        String? fastp_docker
    }
    call fastp_stats {
        input:
            fastq_r1 = fastq_r1,
            fastq_r2 = fastq_r2,
            prefix = prefix,
            fastp_docker = fastp_docker
    }
    output {
        File trimm_r1 = fastp_stats.trimm_r1
        File trimm_r2 = fastp_stats.trimm_r2
        File fastp_html = fastp_stats.fastp_html
        File fastp_json = fastp_stats.fastp_json
    }
}