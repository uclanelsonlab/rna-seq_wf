process samtools_view {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    cpus 12
    tag "Samtools view on $meta"

    input:
    val meta
    path bwa_bam
    val reference

    output:
    path "${meta}_${reference}.view.bam", emit: view_bam
    
    script:
    """
    samtools view -@ $task.cpus -bS -F 2304 -o ${meta}_${reference}.view.bam ${bwa_bam}
    """
}

process samtools_flagstat {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    tag "Samtools flagstat on $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    val meta
    path view_bam
    val reference

    output:
    path "${meta}_${reference}.flagstat.txt", emit: view_bam
    
    script:
    """
    samtools flagstat ${view_bam} > ${meta}_${reference}.flagstat.txt
    """

}

process samtools_index {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    tag "Samtools index on $bam"

    input:
    path reads_gene
    path reads_gene_log
    path final_log
    path sj_tab
    path bam

    output:
    path "${bam}.bai", emit: bam_index
    
    script:
    """
    samtools index ${bam}
    """
}

process samtools_view_sj {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    cpus 12
    tag "Samtools view on $meta for stdout"
    publishDir params.outdir, mode:'symlink'

    input:
    val meta
    path reads_gene
    path reads_gene_log
    path final_log
    path sj_tab
    path star_bam

    output:
    path "tmp_${meta}_bamview.sam", emit: sam_view
    
    script:
    """
    samtools view -@ $task.cpus -h ${star_bam} > tmp_${meta}_bamview.sam
    """
}