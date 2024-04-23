process check_star_reference {
    tag "Check the STAR index to download"

    input:
    tuple val(meta), path(reads)

    output:
    path "star_index", emit: star_index
    val sjdb_overhang
    
    script:
    """
    read_length=\$(zcat ${reads[0]} | head -n 1000 | wc --max-line-length)
    case \$read_length in
        151)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/star_index_151bp_gencode19.tar"
            sjdb_overhang=150
            ;;
        150)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/star2.6.0c_index_150bp_gencode19.tar"
            sjdb_overhang=149
            ;;
        120)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/star_index_120bp_gencode19.tar"
            sjdb_overhang=119
            ;;
        100)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/star_index_100bp_gencode19.tar"
            sjdb_overhang=99
            ;;
        75)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/star_index_75bp_gencode19.tar"
            sjdb_overhang=74
            ;;
        69)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/star_index_69bp_gencode19.tar"
            sjdb_overhang=68
            ;;
        *)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/star_index_100bp_gencode19.tar"
            sjdb_overhang=99
            ;;
    esac

    aws s3 cp \${star_index} .
    mkdir -p star_index/ && tar -xvf \$(basename \${star_index}) -C star_index/
    """
}