process RNASEQC {
    container "gcr.io/broad-cga-aarong-gtex/rnaseqc:latest"
    cpus 40
    tag "Running RNA-SeQC on $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    path gencode_gtf
    tuple val(meta), path(bam)
    tuple val(meta), path(log)
    path(versions)

    output:
    path "*.coverage.tsv",        emit: coverage
    path "*.exon_cv.tsv",         emit: exon_cv
    path "*.exon_reads.gct",      emit: exon_reads
    path "*.gene_fragments.gct",  emit: gene_fragments
    path "*.gene_reads.gct",      emit: gene_reads 
    path "*.gene_tpm.gct",        emit: gene_tpm
    path "*.metrics.tsv",         emit: metrics
    path '*.log',                 emit: log
    path "*versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    rnaseqc ${gencode_gtf} ${bam} --coverage . 2> >(tee ${prefix}.rnaseqc.log >&2)

    cat <<-END_VERSIONS > rnaseqc_versions.yml
    "${task.process}":
        rnaseqc: \$(echo \$(rnaseqc --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}