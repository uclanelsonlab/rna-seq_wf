version 1.0

task samtools_index {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for SAMtools index BAM files"
    }
    input {
        File bam
        String? samtools_docker
    }
    String actual_samtools_docker=select_first([samtools_docker, "quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2"])
    command {
        set -uexo pipefail
        samtools index -@ 12 ~{bam}
    }
    runtime {
        docker: actual_samtools_docker
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
        dx_timeout: "2H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File bai = "${bam}.bai"
    }
}

workflow samtools_index_wf {
    parameter_meta {
        bam: {
            description: "BAM aligned by STAR alignReads",
            extension: ".bam"
        }
        samtools_docker: {
            description: "TAR zip docker image from SAMtools stored on DNAnexus ()",
            extension: ".tar.gz"
        }
    }
    input {
        File bam
        String? samtools_docker
    }
    call samtools_index {
        input:
            bam=bam,
            samtools_docker=samtools_docker
    }
    output {
        File bai = samtools_index.bai
    }
}