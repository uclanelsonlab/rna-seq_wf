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

process download_files { 
    tag "Download ${sample_name} FASTQ and reference files"
    // publishDir params.outdir, mode:'symlink'

    input: 
    val sample_name
    val library
    val fastq_bucket
    val rib_reference_path

    output:
    path "${sample_name}*.fastq.gz", emit: fastq_files
    path "human_rRNA_strict*", emit: rrna_reference_files
    path "human_globinRNA*", emit: globinrna_reference_files

    script: 
    """
    aws s3 cp ${fastq_bucket}/${library}/ . --exclude "*" --recursive --include "${sample_name}*"
    aws s3 cp ${rib_reference_path}/rrna/ . --recursive
    aws s3 cp ${rib_reference_path}/globinrna/ . --recursive
    """
}

workflow { 
    download_files_ch = download_files(params.sample_name, params.library, params.fastq_bucket, params.rib_reference_path)
    download_files_ch.view()
}