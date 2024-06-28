process upload_files {
    tag "Upload all the necessary output files"

    input:
    val library
    val sample_name
    val output_bucket
    path flagstat_rrna
    path flagstat_globinrna
    path reads_gene
    path reads_gene_log
    path final_log
    path sj_tab
    path star_bam
    path sj_tab_gz
    path rare_junctions
    path all_rare_junctions
    path gene_counts
    path gene_counts_short
    path gene_counts_summary
    path rna_cram
    path rna_crai

    script:
    """
    aws s3 cp ${flagstat_rrna} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${flagstat_globinrna} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${reads_gene} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${reads_gene_log} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${final_log} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${sj_tab} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${sj_tab_gz} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${all_rare_junctions} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${rare_junctions} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${gene_counts} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${gene_counts_short} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${gene_counts_summary} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${rna_cram} ${output_bucket}/${sample_name}/${library}/
    aws s3 cp ${rna_crai} ${output_bucket}/${sample_name}/${library}/
    """
}