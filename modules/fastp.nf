process run_fastp {
    container "quay.io/biocontainers/fastp:0.23.3--h5f740d0_0"
    cpus 36
    tag "Fastp on $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz'),  emit: reads
    tuple val(meta), path('*.json'),            emit: json
    tuple val(meta), path('*.html'),            emit: html
    tuple val(meta), path('*.log'),             emit: log
    path "*versions.yml",                       emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    fastp -w $task.cpus -i ${reads[0]} -I ${reads[1]} -o ${prefix}_R1.fastp.fastq.gz -O ${prefix}_R2.fastp.fastq.gz -j ${prefix}.fastp.json -h ${prefix}.fastp.html --detect_adapter_for_pe 2> >(tee ${prefix}.fastp.log >&2)
    
    cat <<-END_VERSIONS > fastp_versions.yml
    "${task.process}":
        fastp: \$(echo \$(fastp --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}