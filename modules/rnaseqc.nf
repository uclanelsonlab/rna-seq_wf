process RNASEQC {
    container "gcr.io/broad-cga-aarong-gtex/rnaseqc:latest"
    cpus 40
    tag "Running RNA-SeQC on $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    val meta
    path gencode_gtf
    path bam
    path versions

    output:
    tuple val(meta), path("*.coverage.tsv"),        emit: coverage
    tuple val(meta), path("*.exon_cv.tsv"),         emit: exon_cv
    tuple val(meta), path("*.exon_reads.gct"),      emit: exon_reads
    tuple val(meta), path("*.gene_fragments.gct"),  emit: gene_fragments
    tuple val(meta), path("*.gene_reads.gct"),      emit:gene_reads 
    tuple val(meta), path("*.gene_tpm.gct"),        emit: gene_tpm
    tuple val(meta), path("*.metrics.tsv"),         emit: metrics
    tuple val(meta), path('*.log'),                 emit: log
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rnaseqc ${gencode_gtf} ${bam} --coverage . 2> >(tee ${prefix}.rnaseqc.log >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaseqc: \$(echo \$(rnaseqc --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}