process bwa_mem {
    container "biocontainers/bwa:v0.7.17_cv1"
    cpus 10
    tag "BWA mem on $reference for $meta"

    input:
    tuple val(meta), path(reads)
    path reference_dir
    val fasta
    val reference

    output:
    path "${meta}_${reference}.bwa.bam", emit: bwa_bam
    
    script:
    """
    bwa mem -t $task.cpus ${reference_dir}/${fasta} ${reads[0]} ${reads[1]} > ${meta}_${reference}.bwa.bam
    """
}