// sudo apt-get install -y libhdf5-dev
process prioritize_splice_junctions {
    tag "Create SJ from BAM for $meta"
    conda 'anaconda::pandas=1.5.3 conda-forge::openpyxl anaconda::pytables'
    publishDir params.outdir, mode:'symlink'

    input:
    val meta
    path sj_tab_gz
    val latest_directory_date
    val latest_udn_id_key
    val sjdblist
    val genemap
    val constraint
    val gencode_cds

    output:
    path "${meta}_rare_junctions_filtered.xlsx", emit: rare_junctions
    path "${meta}_rare_junctions_all.tsv", emit: all_rare_junctions

    script:
    """
    export HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/
    aws s3 cp ${latest_directory_date} nelson_lab_splice_junction_database_hdf5/ --recursive
    aws s3 cp ${latest_udn_id_key} latest_udn_id_key.tsv
    aws s3 cp ${sjdblist} sjdbList.fromGTF.out.tab
    aws s3 cp ${genemap} genemap2_2020-08-02_disease_genes_with_hg19.txt
    aws s3 cp ${constraint} constraint_with_interval_canonical.txt
    aws s3 cp ${gencode_cds} gencode_exons_cds_only.tsv

    for chrom in {1..22} X Y; do
        get_rare_junctions_plus_1_dnanexus.py \
        --proband_bam2sj ${sj_tab_gz} \
        --chromosome \${chrom} \
        --reference_bam2sj_hdf5 nelson_lab_splice_junction_database_hdf5/nelson_lab_sjdb_\${chrom}.h5.gz \
        --udn_id_key latest_udn_id_key.tsv \
        --output_dir . --sjdb sjdbList.fromGTF.out.tab
    done

    tsv_filename=${meta}_rare_junctions_all.tsv
    head -n 1 ${meta}_rare_junctions_chr_22.tsv > \${tsv_filename}
    for chrom in {1..22} X Y; do
        tail -n +2 ${meta}_rare_junctions_chr_\${chrom}.tsv >> \${tsv_filename}
    done

    filtered_tsv_filename=${meta}_rare_junctions_filtered.tsv
    head -n 1 \${tsv_filename} > \${filtered_tsv_filename}
    awk '(\$7 >= 0.4 && \$13 <= 2) || (\$7 >= 0.3 && \$12 <= 2) || (\$7 >= 0.2 && \$11 <= 2) || (\$7 >= 0.1 && \$10 <= 2) {print}' \${tsv_filename} >> \${filtered_tsv_filename}
    
    add_omim_download_to_tsv.py --input \${filtered_tsv_filename} --omim genemap2_2020-08-02_disease_genes_with_hg19.txt --output_dir .
    
    add_constraints_to_tsv.py --input ${meta}_rare_junctions_filtered_with_omim.tsv --output_dir . --constraint constraint_with_interval_canonical.txt
    xlsx_filename=${meta}_rare_junctions_filtered.xlsx

    postprocess.py --input_tsv ${meta}_rare_junctions_filtered_with_omim_with_constraints.tsv --output_filename \${xlsx_filename} --gencode gencode_exons_cds_only.tsv
    """
}