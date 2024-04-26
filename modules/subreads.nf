process gtf_download {
    tag "Download reference GTF file for subread featureCounts"   

    input:
    path gencode_gtf_path

    output:
    path "${meta}.gene_id.exon.ct", emit: gene_counts
    path "${meta}.gene_id.exon.ct.short.txt", emit: gene_counts_short
    path "${meta}.gene_id.exon.ct.summary", emit: gene_counts_summary

    script:
    """
    aws s3 cp ${gencode_gtf_path} .s
    """
}

process subread_featurecounts {
    container "quay.io/biocontainers/subread:2.0.6--he4a0461_0"
    cpus 32
    tag "Subreads featureCounts on $meta"   

    input:
    val meta
    path gencode_gtf
    path rna_bam

    output:
    path "${meta}.gene_id.exon.ct", emit: gene_counts
    path "${meta}.gene_id.exon.ct.short.txt", emit: gene_counts_short
    path "${meta}.gene_id.exon.ct.summary", emit: gene_counts_summary

    script:
    """
    featureCounts -T 4 -t exon -g gene_id -a ${gencode_gtf} -o ${meta}.gene_id.exon.ct -p -C --primary ${rna_bam}
    awk -F $' ' 'BEGIN {OFS=FS} { print \$1, \$7 }' ${meta}.gene_id.exon.ct > ${meta}.gene_id.exon.ct.short.txt
    """
}