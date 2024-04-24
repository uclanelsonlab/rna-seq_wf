#!/usr/bin/env nextflow

params.sample_name = 'NB-8204-M-muscle'
params.library = "SN_7RNA_S-24-0479_XA044"
params.fastq_bucket = "s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq"
params.rib_reference_path = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference"
params.outdir = "results"

log.info """\
    R N A - S E Q _ W F   P I P E L I N E
    ===================================
    sample_name         : ${params.sample_name}
    library             : ${params.library}
    fastq_bucket        : ${params.fastq_bucket}
    rib_reference_path  : ${params.rib_reference_path}
    outdir              : ${params.outdir}
    """
    .stripIndent(true)

include { download_fastqs; download_rrna; download_globinrna   } from './modules/download_files.nf'
include { run_fastp } from './modules/fastp.nf'
include { filter_fastq } from './modules/filters.nf'
include { bwa_mem as bwa_mem_rrna } from './modules/bwa.nf'
include { bwa_mem as bwa_mem_globinrna } from './modules/bwa.nf'
include { samtools_view as samtools_view_rrna; samtools_flagstat as samtools_flagstat_rrna; samtools_index } from './modules/samtools.nf'
include { samtools_view as samtools_view_globinrna; samtools_flagstat as samtools_flagstat_globinrna } from './modules/samtools.nf'
include { check_star_reference; star_alignreads } from './modules/star.nf'

workflow {
    download_fastqs_ch = download_fastqs(params.sample_name, params.library, params.fastq_bucket)
    download_rrna_ch = download_rrna(params.rib_reference_path)
    download_globinrna_ch = download_globinrna(params.rib_reference_path)

    // contamination check
    fastp_ch = run_fastp(download_fastqs_ch)
    filtered_fastq_ch = filter_fastq(fastp_ch)
    rrna_bwa_ch = bwa_mem_rrna(filtered_fastq_ch, download_rrna_ch, "human_rRNA_strict.fasta", "rrna")
    globinrna_bwa_ch = bwa_mem_globinrna(filtered_fastq_ch, download_globinrna_ch, "human_globinRNA.fa", "globinrna")
    rrna_samtools_view_ch = samtools_view_rrna(params.sample_name, rrna_bwa_ch, "rrna")
    globinrna_samtools_view_ch = samtools_view_globinrna(params.sample_name, globinrna_bwa_ch, "globinrna")
    rrna_samtools_flagstat_ch = samtools_flagstat_rrna(params.sample_name, rrna_samtools_view_ch, "rrna")
    globinrna_samtools_flagstat_ch = samtools_flagstat_globinrna(params.sample_name, globinrna_samtools_view_ch, "globinrna")

    // STAR alignment
    star_index_ref_ch = check_star_reference(download_fastqs_ch)
    star_alignreads_ch = star_alignreads(params.sample_name, star_index_ref_ch, fastp_ch)
}