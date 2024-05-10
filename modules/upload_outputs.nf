process upload_files {
    tag "Upload all the necessary output files"

    input:
    val library
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
    aws s3 cp ${flagstat_rrna} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/star_output/${library}/${flagstat_rrna}
    aws s3 cp ${flagstat_globinrna} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/star_output/${library}/${flagstat_globinrna}
    aws s3 cp ${reads_gene} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/star_output/${library}/${reads_gene}
    aws s3 cp ${reads_gene_log} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/star_output/${library}/${reads_gene_log}
    aws s3 cp ${final_log} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/star_output/${library}/${final_log}
    aws s3 cp ${sj_tab} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/star_output/${library}/${sj_tab}
    aws s3 cp ${sj_tab_gz} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/bam2sj/${library}/${sj_tab_gz}
    aws s3 cp ${all_rare_junctions} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/rare_junctions/${library}/${all_rare_junctions}
    aws s3 cp ${rare_junctions} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/rare_junctions/${library}/${rare_junctions}
    aws s3 cp ${gene_counts} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/featurecounts_pc_unstranded/${library}/${gene_counts}
    aws s3 cp ${gene_counts_short} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/featurecounts_pc_unstranded/${library}/${gene_counts_short}
    aws s3 cp ${gene_counts_summary} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/featurecounts_pc_unstranded/${library}/${gene_counts_summary}
    aws s3 cp ${rna_cram} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/star_output/${library}/${rna_cram}
    aws s3 cp ${rna_crai} s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/star_output/${library}/${rna_crai}
    """
}