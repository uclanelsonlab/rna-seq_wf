version 1.0

task qualimap_rnaseq {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Qualimap, RNA-seq QC reports quality control metrics and bias estimations which are specific for whole transcriptome sequencing"
    }
    input {
        File bam_file
        File bam_index
        File gtf
        String? mem
        String? qualimap_docker
    }

    Int actual_mem=select_first([mem, "20G"])
    String actual_qualimap_docker=select_first([qualimap_docker, "quay.io/biocontainers/qualimap:2.2.2d--hdfd78af_2"])
    command {
        set -uexo pipefail
        ls ~{bam_index}
        qualimap rnaseq \
            --java-mem-size=~{actual_mem} \
            -outdir qualimap/ \
            -a proportional \
            -bam ~{bam_file} \
            -p strand-specific-reverse \
            -gtf ~[gtf]
    }
    runtime {
        docker: actual_qualimap_docker
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
        File qualimapReport = "qualimap/qualimapReport.html"
        File qualimap_results = "qualimap/rnaseq_qc_results.txt"
    }
}

workflow qualimap_rnaseq_wf {
    parameter_meta {
        bam_file: {
            description: "BAM file",
            extension: ".bam"
        }
        bam_index: {
            description: "BAM index file",
            extension: ".bai"
        }
        qualimap_docker: {
            description: "TAR zip docker image from Qualimap stored on DNAnexus (file-GVzvBQ802k8k4zPYKqpfz1vZ)",
            extension: ".tar.gz"
        }
        gtf: {
            description: "GTF annotation for reference genome"
        }
        mem: {
            description: "Memory string to be used, default: 20G"
        }
    }
    input {
        File bam_file
        File bam_index
        File gtf
        String? mem
        String? qualimap_docker
    }
    call qualimap_rnaseq {
        input:
            bam_file = bam_file,
            bam_index = bam_index,
            gtf = gtf,
            mem = mem,
            qualimap_docker = qualimap_docker
    }
    output {
        File qualimapReport = qualimap_rnaseq.qualimapReport
        File qualimap_results = qualimap_rnaseq.qualimap_results
    }
}