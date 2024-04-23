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

process download_rrna {
    tag "Download rRNA files"

    input:
    val rib_reference_path

    output:
    path "rrna_reference", emit: rrna_reference_dir

    script:
    """
    mkdir rrna_reference
    aws s3 cp ${rib_reference_path}/rrna/ rrna_reference/ --recursive
    """
}

process download_globinrna {
    tag "Download globinRNA files"

    input:
    val rib_reference_path

    output:
    path "globinrna_reference", emit: globinrna_reference_dir

    script:
    """
    mkdir globinrna_reference
    aws s3 cp ${rib_reference_path}/globinrna/ globinrna_reference/ --recursive
    """
}