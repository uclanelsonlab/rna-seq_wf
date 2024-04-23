process filter_fastq {
    tag "Filter $meta FASTQ files"
    // publishDir params.outdir, mode:'symlink'

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(json)
    tuple val(meta), path(html)
    tuple val(meta), path(log)

    output:
    tuple val(meta), path(reads), emit: reads

    script:
    """
    zcat ${reads[0]} | head -n 4000000 | gzip > ${meta}_R1.4kreads.fastq.gz
    zcat ${reads[1]} | head -n 4000000 | gzip > ${meta}_R2.4kreads.fastq.gz
    """
}