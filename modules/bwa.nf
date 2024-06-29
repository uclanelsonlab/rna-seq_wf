process bwa_mem {
    container "biocontainers/bwa:v0.7.17_cv1"
    cpus 12
    tag "BWA mem on $reference for $meta"

    input:
    tuple val(meta), path(reads)
    path reference_dir
    val fasta
    val reference

    output:
    tuple val(meta), path("*.bwa.bam"), emit: bwa_bam
    tuple val(meta), path('*.log'),     emit: log
    path "versions.yml",                emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bwa mem -t $task.cpus ${reference_dir}/${fasta} ${reads[0]} ${reads[1]} > ${prefix}_${reference}.bwa.bam 2> >(tee ${prefix}.bwa.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}