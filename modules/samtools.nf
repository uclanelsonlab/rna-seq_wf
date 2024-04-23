process samtools_view {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    cpus 10
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