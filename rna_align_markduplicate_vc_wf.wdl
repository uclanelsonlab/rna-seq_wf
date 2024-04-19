version 1.0

import "tasks/fastp.wdl" as fastp_task
import "tasks/star_align_reads.wdl" as star_align_reads_task
import "tasks/samtools_index.wdl" as samtools_index_task
import "tasks/picard_markduplicates.wdl" as picard_markduplicates_task
import "tasks/picard.wdl" as picard_task
import "tasks/qualimap.wdl" as qualimap_task
import "tasks/multiqc.wdl" as multiqc_task
import "tasks/gatk.wdl" as gatk_task

workflow rna_align_markduplicate_vc_wf {
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
        fastp_docker: {
            description: "TAR zip docker image from fastp stored on DNAnexus (file-)",
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
        File fasta
        File fasta_fai
        File fasta_dict
        # fastp/star_align_reads
        String prefix
        File fastq_r1
        File fastq_r2
        File star_index
        String? fastp_docker
        String? star_docker
        String? gatk_docker
        String? samtools_docker
        String? multiqc_image
        # STAR options
        Int? runThreadN
        Int? sjdbOverhang
        Int? limitBAMsortRAM
        Int? outBAMsortingThreadN
        # Qualimap
        File gtf
        String? java_mem
        String? mem
        String? qualimap_docker
        # gatk
        File dbsnp138
        File known_indels
        File indels_1kG
        File af_only_gnomad
        File small_exac_common_3
    }

    call fastp_task.fastp_stats as fastp_stats {
        input:
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            prefix=prefix,
            fastp_docker=fastp_docker
    }

    call star_align_reads_task.align_reads as align_reads {
        input:
            fastq_r1=fastp_stats.trimm_r1,
            fastq_r2=fastp_stats.trimm_r2,
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

    call qualimap_task.qualimap_rnaseq as qualimap_rnaseq {
        input:
            bam_file=picard_marked_dup.marked_bam,
            bam_index=samtools_filter_2.bai,
            gtf=gtf,
            mem=mem,
            qualimap_docker=qualimap_docker
    }
    call picard_task.picard_CollectMultipleMetrics as picard_CollectMultipleMetrics {
        input:
            fasta = fasta,
            fasta_fai = fasta_fai,
            fasta_dict = fasta_dict,
            bam_file = picard_marked_dup.marked_bam,
            bam_index = samtools_filter_2.bai,
            prefix = prefix,
            java_mem = java_mem,
            picard_docker = gatk_docker
    }

    call gatk_task.gatk_SplitNCigarReads as gatk_SplitNCigarReads {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            bam=picard_marked_dup.marked_bam,
            bai=samtools_filter_2.bai,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_task.gatk_BaseRecalibrator {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            bam=gatk_SplitNCigarReads.bam_split_cigar,
            bai=gatk_SplitNCigarReads.bai_split_cigar,
            dbsnp138=dbsnp138,
            known_indels=known_indels,
            indels_1kG=indels_1kG,
            af_only_gnomad=af_only_gnomad,
            small_exac_common_3=small_exac_common_3,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_task.gatk_ApplyBQSR as gatk_ApplyBQSR {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            bam=gatk_SplitNCigarReads.bam_split_cigar,
            bai=gatk_SplitNCigarReads.bai_split_cigar,
            recal_data=gatk_BaseRecalibrator.recal_data,
            prefix=prefix,
            gatk_docker=gatk_docker
    }    

    call gatk_task.gatk_AnalyzeCovariates as gatk_AnalyzeCovariates {
        input:
            recal_data=gatk_BaseRecalibrator.recal_data,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_task.gatk_HaplotypeCaller as gatk_HaplotypeCaller {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            bam=gatk_ApplyBQSR.bqsr_bam,
            bai=gatk_ApplyBQSR.bqsr_bai,
            dbsnp138=dbsnp138,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_task.gatk_GenotypeGVCFs as gatk_GenotypeGVCFs {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            hc_gvcf=gatk_HaplotypeCaller.hc_gvcf,
            hc_gvcf_index=gatk_HaplotypeCaller.hc_gvcf_index,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_task.gatk_VariantFiltration as gatk_VariantFiltration {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            hc_vcf=gatk_GenotypeGVCFs.hc_vcf,
            hc_vcf_index=gatk_GenotypeGVCFs.hc_vcf_index,
            prefix=prefix,
            gatk_docker=gatk_docker
    }
    call multiqc_task.multiqc_array as multiqc {
        input:
            stats_files = [
                fastp_stats.fastp_html, 
                fastp_stats.fastp_json, 
                qualimap_rnaseq.qualimapReport,
                qualimap_rnaseq.qualimap_results,
                gatk_AnalyzeCovariates.analyze_pdf,
                picard_CollectMultipleMetrics.alignment_summary_metrics,
                picard_CollectMultipleMetrics.base_distribution_by_cycle,
                picard_CollectMultipleMetrics.base_distribution_by_cycle_metrics,
                picard_CollectMultipleMetrics.insert_size_histogram,
                picard_CollectMultipleMetrics.insert_size_metrics,
                picard_CollectMultipleMetrics.quality_by_cycle,
                picard_CollectMultipleMetrics.quality_by_cycle_metrics,
                picard_CollectMultipleMetrics.quality_distribution,
                picard_CollectMultipleMetrics.quality_distribution_metrics,
                picard_CollectMultipleMetrics.read_length_histogram
            ],
            prefix = prefix,
            multiqc_image = multiqc_image
    }
    output {
        File star_bam = align_reads.star_bam
        File sjdb_txt = align_reads.sjdb_txt
        File sjdb_tab = align_reads.sjdb_tab
        File read_counts = align_reads.read_counts
        File junctions = align_reads.junctions
        File junctions_pass1 = align_reads.junctions_pass1
        File junctions_pass1_log = align_reads.junctions_pass1_log
        Array[File] logs = align_reads.logs
        File bai = samtools_filter_2.bai
        File qualimapReport = qualimap_rnaseq.qualimapReport
        File qualimap_results = qualimap_rnaseq.qualimap_results 
        File bqsr_bam = gatk_ApplyBQSR.bqsr_bam
        File bqsr_bai = gatk_ApplyBQSR.bqsr_bai
        File analyze_pdf = gatk_AnalyzeCovariates.analyze_pdf
        File hc_gvcf = gatk_HaplotypeCaller.hc_gvcf
        File hc_gvcf_index = gatk_HaplotypeCaller.hc_gvcf_index
        File hc_vcf = gatk_GenotypeGVCFs.hc_vcf
        File hc_vcf_index = gatk_GenotypeGVCFs.hc_vcf_index
        File filtered_vcf = gatk_VariantFiltration.filtered_vcf
        File filtered_vcf_index = gatk_VariantFiltration.filtered_vcf_index
        File multiqc_html = multiqc.multiqc_html
        File multiqc_data = multiqc.multiqc_data
    }
}