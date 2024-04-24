process check_star_reference {
    tag "Check the STAR index to download"

    input:
    tuple val(meta), path(reads)

    output:
    path "star_index", emit: star_index
    env sjdb_overhang
    
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

process star_alignreads {
    container "gvcn/star:2.6.0c"
    cpus 36
    tag "STAR alignReadson $meta"   

    input:
    val meta
    val reference
    val sjdb_overhang
    tuple val(meta), path(reads)
    tuple val(meta), path(json)
    tuple val(meta), path(html)
    tuple val(meta), path(log)

    output:
    path "./star_output/${meta}.ReadsPerGene.out.tab.gz", emit: reads_gene
    path "./star_output/${meta}.ReadsPerGene.log.out", emit: reads_gene_log
    path "./star_output/${meta}.Log.final.out", emit: final_log
    path "./star_output/${meta}.SJ.out.tab.gz", emit: sj_tab
    path "./star_output/${meta}.Aligned.sortedByCoord.out.bam", emit: star_bam

    script:
    """
    sort_mem=\$(head -n1 /proc/meminfo | awk '{print int(\$2*1000*0.90)}')
    STAR --runMode alignReads --runThreadN $task.cpus --genomeDir ${reference} --twopassMode Basic --sjdbOverhang ${sjdb_overhang} --readFilesIn ${reads[0]} ${reads[1]} --readFilesCommand zcat --outFileNamePrefix ${meta} --alignSoftClipAtReferenceEnds Yes --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outBAMcompression -1 --outSAMunmapped Within --genomeLoad NoSharedMemory --limitBAMsortRAM \$sort_mem --outBAMsortingThreadN $task.cpus --outSAMattrRGline ID:rg1 SM:${meta} PL:Illumina LB:${meta}
    """
}