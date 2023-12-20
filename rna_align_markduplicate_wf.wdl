version 1.0

import "tasks/star_align_reads.wdl" as star_align_reads_task
import "tasks/samtools_index.wdl" as samtools_index_task
import "tasks/picard_markduplicates.wdl" as picard_markduplicates_task

workflow rna_align_markduplicate_wf {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "# Workflow for RNA-Seq alignment"
    }
    parameter_meta {
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        fastq_r1: {
            description: "Sample R1 FASTQ file",
            extension: ".R1.fastq.gz"
        }
        fastq_r2: {
            description: "Sample R2 FASTQ file",
            extension: ".R2.fastq.gz"
        }
        star_index: {
            description: "TAR zip reference files from Gencode",
            extension: ".tar.gz"
        }
        star_docker: {
            description: "TAR zip docker image from STAR stored on DNAnexus (file-GVjJgfQ02k8bYxZ9z70g869k)",
            extension: ".tar.gz"
        }
        gatk_docker: {
            description: "GATK docker image (file-GVgyqvQ02k8XGzgY2V8BVVBG)"
        }
        samtools_docker: {
            description: "TAR zip docker image from SAMtools stored on DNAnexus ()",
            extension: ".tar.gz"
        }
    }
    input {
        # star_align_reads
        String prefix
        File fastq_r1
        File fastq_r2
        File star_index
        String? star_docker
        String? gatk_docker
        # STAR options
        Int? runThreadN
        Int? sjdbOverhang
        Int? limitBAMsortRAM
        Int? outBAMsortingThreadN

        # samtools_vw_filter
        String? samtools_docker
    }

    call star_align_reads_task.align_reads as align_reads {
        input:
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            prefix=prefix,
            star_index=star_index,
            star_docker=star_docker
    }
    call samtools_index_task.samtools_index as samtools_filter_1 {
        input:
            bam=align_reads.star_bam,
            samtools_docker=samtools_docker
    }
    call picard_markduplicates_task.picard_marked_dup as picard_marked_dup {
        input:
            bam=align_reads.star_bam,
            bai=samtools_filter_1.bai,
            prefix=prefix,
            picard_docker=gatk_docker
    }
    call samtools_index_task.samtools_index as samtools_filter_2 {
        input:
            bam=picard_marked_dup.marked_bam,
            samtools_docker=samtools_docker
    }

    output {
        # align_reads
        File star_bam = align_reads.star_bam
        File sjdb_txt = align_reads.sjdb_txt
        File sjdb_tab = align_reads.sjdb_tab
        File read_counts = align_reads.read_counts
        File junctions = align_reads.junctions
        File junctions_pass1 = align_reads.junctions_pass1
        File junctions_pass1_log = align_reads.junctions_pass1_log
        Array[File] logs = align_reads.logs
        # samtools_vw_filter
        File bam = picard_marked_dup.marked_bam
        File bai = samtools_filter_2.bai
    }
}