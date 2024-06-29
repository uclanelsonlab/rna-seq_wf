
process download_gencode {
    tag "Download reference GTF file for subread featureCounts"

    input:
    val gencode_gtf_path

    output:
    path("*.gtf"), emit: gencode_gtf

    script:
    """
    aws s3 cp ${gencode_gtf_path} .
    """
}

process subread_featurecounts {
    tag "Generate counts by gene using featureCounts"
    container "quay.io/biocontainers/subread:2.0.6--he4a0461_0"
    cpus 12
    publishDir params.outdir, mode:'symlink'

    input:
    path gencode_pc
    tuple val(meta), path(bam)
    tuple val(meta), path(log)
    path(versions)

    output:
    path "*.gene_id.exon.ct",             emit: gene_counts
    path "*.gene_id.exon.ct.short.txt",   emit: gene_counts_short
    path "*.gene_id.exon.ct.summary",     emit: gene_counts_summary
    path "*.log",                         emit: log
    path "*versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    featureCounts -T $task.cpus -t exon -g gene_id -a ${gencode_pc} -o ${prefix}.gene_id.exon.ct -p -C --primary ${bam} 2> >(tee ${prefix}.featureCounts.log >&2)
    awk -F \$' ' 'BEGIN {OFS=FS} { print \$1, \$7 }' ${prefix}.gene_id.exon.ct > ${prefix}.gene_id.exon.ct.short.txt
    
    cat <<-END_VERSIONS > featureCounts_versions.yml
    "${task.process}":
        featureCounts: \$(echo \$(featureCounts -v 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}