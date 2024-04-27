process bam2sj {
    tag "Create SJ from BAM for $meta"
    container "gvcn/bam2sj:v0.0.2"
    cache false

    input:
    path meta
    path sam_view

    output:
    path "${meta}.bam2SJ.out.tab.gz", emit: sj_tab_gz

    script:
    """
    perl bam2sj.pl ${sam_view} | sort -k1,1 -k2,2n > ${meta}.bam2SJ.out.tab
    gzip ${meta}.bam2SJ.out.tab
    """
}