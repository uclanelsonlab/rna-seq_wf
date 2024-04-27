process bam2sj {
    tag "Create SJ from BAM for $meta"

    input:
    val meta
    path sam_view

    output:
    path "${meta}.bam2SJ.out.tab.gz", emit: sj_tab_gz

    script:
    """
    bam2sj.pl ${sam_view} > ${meta}.bam2SJ.out.tab
    gzip ${meta}.bam2SJ.out.tab
    """
}