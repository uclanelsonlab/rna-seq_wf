process upload_files {
    tag "Upload all the necessary output files"

    input:
    val library
    val sample_name
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

    script:
    """
    aws s3 cp ${flagstat_rrna} ${output_bucket}/${sample_name}/${library}/${flagstat_rrna}
    aws s3 cp ${flagstat_rrna_v} ${output_bucket}/${sample_name}/${library}/${flagstat_rrna_v}
    aws s3 cp ${flagstat_globinrna} ${output_bucket}/${sample_name}/${library}/${flagstat_globinrna}
    aws s3 cp ${flagstat_globinrna_v} ${output_bucket}/${sample_name}/${library}/${flagstat_globinrna_v}
    aws s3 cp ${reads_gene} ${output_bucket}/${sample_name}/${library}/${reads_gene}
    aws s3 cp ${reads_gene_log} ${output_bucket}/${sample_name}/${library}/${reads_gene_log}
    aws s3 cp ${final_log} ${output_bucket}/${sample_name}/${library}/${final_log}
    aws s3 cp ${sj_tab} ${output_bucket}/${sample_name}/${library}/${sj_tab}
    aws s3 cp ${star_log} ${output_bucket}/${sample_name}/${library}/${star_log}
    """
}