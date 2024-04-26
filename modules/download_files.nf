process download_fastqs {
    tag "Download ${meta} FASTQ files"

    input:
    val meta
    val library
    val fastq_bucket

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads

    script:
    """
    aws s3 cp ${fastq_bucket}/${library}/ . --exclude "*" --recursive --include "${meta}*"
    """
}

process download_rna_ref {
    tag "Download rna reference files"

    input:
    val rna_reference_path
    val type

    output:
    path "${type}_reference", emit: rrna_reference_dir

    script:
    """
    mkdir ${type}_reference
    aws s3 cp ${rna_reference_path}/${type}/ ${type}_reference/ --recursive
    """
}