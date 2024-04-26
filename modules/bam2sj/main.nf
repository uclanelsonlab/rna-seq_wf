process samtools_view_sj {
    conda "conda-forge::perl bioconda::samtools"
    cpus 12
    tag "Samtools view on $meta for stdout"
    publishDir params.outdir, mode:'symlink'
    cache false

    input:
    val meta
    path reads_gene
    path reads_gene_log
    path final_log
    path sj_tab
    path star_bam

    output:
    path "${meta}.bam2SJ.out.tab.gz", emit: sj_tab_gz
    
    script:
    """
    samtools view -@ $task.cpus -h ${star_bam} | perl bam2sj.pl | sort -k1,1 -k2,2n > ${meta}.bam2SJ.out.tab
    gzip ${meta}.bam2SJ.out.tab
    """
}


workflow bar {
    take:
        value

    main:
        value | samtools_view_sj

    emit:
        samtools_view_sj.out
}