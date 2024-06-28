process SAMBAMBA_MARKDUP {
    container "quay.io/biocontainers/sambamba:1.0.1--h6f6fda4_1"
    cpus 40
    tag "Running sambamba markdup on $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    val meta
    path reads_gene
    path reads_gene_log
    path final_log
    path sj_tab
    path bam

    output:
    path "${meta}.markdup.bam"      , emit: marked_bam
    path "${meta}.markdup.bai"      , emit: marked_bai, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sambamba markdup -t $task.cpus --tmpdir ./ $bam ${meta}.markdup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}