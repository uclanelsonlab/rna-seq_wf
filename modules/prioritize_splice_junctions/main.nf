// sudo apt-get install -y libhdf5-dev
process prioritize_splice_junctions {
    tag "Create SJ from BAM for $meta"
    conda 'anaconda::pandas=1.5.3 conda-forge::openpyxl anaconda::pytables'

    input:
    path sj_tab_gz
    val latest_directory_date
    val latest_udn_id_key
    path sjdblist

    output:
    path "${meta}.bam2SJ.out.tab.gz", emit: sj_tab_gz

    script:
    """
    export HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/
    aws s3 cp ${latest_directory_date} nelson_lab_splice_junction_database_hdf5/ --recursive
    aws s3 cp ${latest_udn_id_key} latest_udn_id_key.tsv

    for chrom in {1..22} X Y; do
        get_rare_junctions_plus_1_dnanexus.py \
        --proband_bam2sj ${sj_tab_gz} \
        --chromosome \${chrom} \
        --reference_bam2sj_hdf5 nelson_lab_splice_junction_database_hdf5/nelson_lab_sjdb_\${chrom}.h5.gz \
        --udn_id_key latest_udn_id_key.tsv \
        --output_dir . --sjdb ${sjdblist}
    done
    """
}