process upload_files {
    tag "Upload all the necessary output files"

    input:
    val library
    val output_bucket
    path flagstat_rrna
    path flagstat_rrna_v
    path flagstat_globinrna 
    path flagstat_globinrna_v
    tuple val(meta), path(reads_gene) 
    tuple val(meta), path(reads_gene_log)
    tuple val(meta), path(final_log)
    tuple val(meta), path(sj_tab)
    tuple val(meta), path(star_bam)
    tuple val(meta), path(star_log)
    path star_versions
    path gene_counts
    path gene_counts_short
    path gene_counts_summary
    path subread_log
    path subread_versions
    path qc_coverage
    path qc_exon_cv
    path qc_exon_reads
    path qc_gene_fragments
    path qc_gene_reads
    path qc_gene_tpm
    path qc_metrics
    path qc_log
    path qc_versions
    tuple val(meta), path(cram)
    tuple val(meta), path(crai)
    path cram_log
    path cram_versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    aws s3 cp ${flagstat_rrna} ${output_bucket}/${prefix}/${library}/qc/${flagstat_rrna}
    aws s3 cp ${flagstat_rrna_v} ${output_bucket}/${prefix}/${library}/qc/${flagstat_rrna_v}
    aws s3 cp ${flagstat_globinrna} ${output_bucket}/${prefix}/${library}/qc/${flagstat_globinrna}
    aws s3 cp ${flagstat_globinrna_v} ${output_bucket}/${prefix}/${library}/qc/${flagstat_globinrna_v}
    aws s3 cp ${reads_gene} ${output_bucket}/${prefix}/${library}/alignment/${reads_gene}
    aws s3 cp ${reads_gene_log} ${output_bucket}/${prefix}/${library}/alignment/${reads_gene_log}
    aws s3 cp ${final_log} ${output_bucket}/${prefix}/${library}/alignment/${final_log}
    aws s3 cp ${sj_tab} ${output_bucket}/${prefix}/${library}/alignment/${sj_tab}
    aws s3 cp ${star_log} ${output_bucket}/${prefix}/${library}/alignment/${star_log}
    aws s3 cp ${star_versions} ${output_bucket}/${prefix}/${library}/alignment/${star_versions}
    aws s3 cp ${gene_counts} ${output_bucket}/${prefix}/${library}/counts/${gene_counts}
    aws s3 cp ${gene_counts_short} ${output_bucket}/${prefix}/${library}/counts/${gene_counts_short}
    aws s3 cp ${gene_counts_summary} ${output_bucket}/${prefix}/${library}/counts/${gene_counts_summary}
    aws s3 cp ${subread_log} ${output_bucket}/${prefix}/${library}/counts/${subread_log}
    aws s3 cp ${subread_versions} ${output_bucket}/${prefix}/${library}/counts/${subread_versions}
    aws s3 cp ${qc_coverage} ${output_bucket}/${prefix}/${library}/qc/${qc_coverage}
    aws s3 cp ${qc_exon_cv} ${output_bucket}/${prefix}/${library}/qc/${qc_exon_cv}
    aws s3 cp ${qc_exon_reads} ${output_bucket}/${prefix}/${library}/qc/${qc_exon_reads}
    aws s3 cp ${qc_gene_fragments} ${output_bucket}/${prefix}/${library}/qc/${qc_gene_fragments}
    aws s3 cp ${qc_gene_reads} ${output_bucket}/${prefix}/${library}/qc/${qc_gene_reads}
    aws s3 cp ${qc_gene_tpm} ${output_bucket}/${prefix}/${library}/qc/${qc_gene_tpm}
    aws s3 cp ${qc_metrics} ${output_bucket}/${prefix}/${library}/qc/${qc_metrics}
    aws s3 cp ${qc_log} ${output_bucket}/${prefix}/${library}/qc/${qc_log}
    aws s3 cp ${qc_versions} ${output_bucket}/${prefix}/${library}/qc/${qc_versions}
    aws s3 cp ${cram} ${output_bucket}/${prefix}/${library}/alignment/${cram}
    aws s3 cp ${crai} ${output_bucket}/${prefix}/${library}/alignment/${crai}
    aws s3 cp ${cram_log} ${output_bucket}/${prefix}/${library}/alignment/${cram_log}
    aws s3 cp ${cram_versions} ${output_bucket}/${prefix}/${library}/alignment/${cram_versions}
    """
}