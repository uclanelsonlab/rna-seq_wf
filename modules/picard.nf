process run_markdup {
    container "broadinstitute/picard:2.27.5"
    cpus 36
    tag "Running Picard MarkDuplicates on $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    val meta
    path reads_gene
    path reads_gene_log
    path final_log
    path sj_tab
    path bam

    output:
    path "${meta}.marked_duplicates.bam", emit: marked_bam
    path "${meta}.marked_dup_metrics.txt", emit: dup_metrics
    
    script:
    """
    java "-Xmx100g" -jar /usr/picard/picard.jar MarkDuplicates -I ${bam} -O ${meta}.marked_duplicates.bam -M ${meta}.marked_dup_metrics.txt
    """
}