version 1.0

import "samtools_index.wdl" as samtools_index_wf

task picard_marked_dup {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK MarkDuplicates"
    }
    input {
        File bam
        File bai
        String prefix
        String? java_mem
        String? picard_docker
    }
    String actual_java_mem=select_first([java_mem, "-Xmx20g"])
    String actual_picard_docker=select_first([picard_docker, "broadinstitute/picard:2.27.5"])
    command {
        set -uexo pipefail
        java ~{actual_java_mem} -jar /usr/picard/picard.jar MarkDuplicates \
            -I ~{bam} \
            -O ~{prefix}.marked_duplicates.bam \
            -M ~{prefix}.marked_dup_metrics.txt
    }
    runtime {
        docker: actual_picard_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "10H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File marked_bam = "${prefix}.marked_duplicates.bam"
        File metrics_bam = "${prefix}.marked_dup_metrics.txt"
    }
}

workflow picard_marked_dup_wf {
    parameter_meta {
        bam: {
            description: "BAM aligned by SAMtools vW filtering",
            extension: ".bam"
        }
        bai: {
            description: "Index for BAM aligned by SAMtools vW filtering",
            extension: ".bai"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        picard_docker: {
            description: "Picard docker image (file-)"
        }
    }
    input {
        File bam
        File bai
        String prefix
        String? picard_docker
        String? samtools_docker
    }
    call picard_marked_dup {
        input: 
            bam=bam,
            bai=bai,
            prefix=prefix,
            picard_docker=picard_docker
    }
    call samtools_index_wf.samtools_index as samtools_index {
        input:
            bam=picard_marked_dup.marked_bam,
            samtools_docker=samtools_docker
    }
    output {
        File marked_bam = picard_marked_dup.marked_bam
        File metrics_bam = picard_marked_dup.metrics_bam
        File marked_bai = samtools_index.bai
    }
}
