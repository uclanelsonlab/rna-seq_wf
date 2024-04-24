process run_fastp {
    container "quay.io/biocontainers/fastp:0.23.3--h5f740d0_0"
    cpus 36
    tag "Fastp on $meta"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    
    script:
    """
    fastp -w $task.cpus -i ${reads[0]} -I ${reads[1]} -o ${meta}_R1.fastp.fastq.gz -O ${meta}_R2.fastp.fastq.gz -j ${meta}.fastp.json -h ${meta}.fastp.html --detect_adapter_for_pe 2> >(tee ${meta}.fastp.log >&2)
    """
}