process subread_featurecounts {
    tag "Generate counts by gene using featureCounts"
    container "quay.io/biocontainers/subread:2.0.6--he4a0461_0"
    cpus 12
    publishDir params.outdir, mode:'symlink'

    input:
    val meta
    path gencode_pc
    path bam
    path dup_metrics

    output:
    path "${meta}.gene_id.exon.ct", emit: gene_counts
    path "${meta}.gene_id.exon.ct.short.txt", emit: gene_counts_short
    path "${meta}.gene_id.exon.ct.summary", emit: gene_counts_summary

    script:
    """
    featureCounts -T $task.cpus -t exon -g gene_id -a ${gencode_pc} -o ${meta}.gene_id.exon.ct -p -C --primary ${bam}
    awk -F \$' ' 'BEGIN {OFS=FS} { print \$1, \$7 }' ${meta}.gene_id.exon.ct > ${meta}.gene_id.exon.ct.short.txt
    """
}

process download_gencode {
    tag "Download reference GTF file for subread featureCounts"

    input:
    val gencode_gtf_path

    output:
    path "gencode.protein_coding.gtf", emit: gencode_pc

    script:
    """
    aws s3 cp ${gencode_gtf_path} gencode.protein_coding.gtf
    """
}