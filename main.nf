nextflow.enable.dsl = 2

params.sample_name = 'NB-8204-M-muscle'
params.library = "SN_7RNA_S-24-0479_XA044"
params.fastq_bucket = "s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq"
params.rib_reference_path = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference"
params.gencode_gtf_path = "s3://ucla-rare-diseases/UCLA-UDN/alden_test/resources/gencode.v19.protein_coding.gtf"
params.outdir = "results"
// params for prioritize splice junctions
params.constraint = "s3://ucla-rare-diseases/UCLA-UDN/assets/prioritize_splice_junctions/constraint_with_interval_canonical.txt"
params.gencode_cds = "s3://ucla-rare-diseases/UCLA-UDN/assets/prioritize_splice_junctions/gencode_exons_cds_only.tsv"
params.genemap = "s3://ucla-rare-diseases/UCLA-UDN/assets/prioritize_splice_junctions/genemap2_2020-08-02_disease_genes_with_hg19.txt"
params.gnomad = "s3://ucla-rare-diseases/UCLA-UDN/assets/prioritize_splice_junctions/gnomad.v2.1.1.lof_metrics.by_gene.txt"
params.sjdblist = "s3://ucla-rare-diseases/UCLA-UDN/assets/prioritize_splice_junctions/sjdbList.fromGTF.out.tab"
params.latest_udn_id_key="s3://ucla-rare-diseases/UCLA-UDN/assets/pedigrees/udn_id_to_proband_id_2020-02-29.tsv" // latest_udn_id_key.tsv
params.latest_directory_date="s3://ucla-rare-diseases/UCLA-UDN/assets/nelson_lab_splice_junction_database_hdf5/2023-09-08/"
params.human_fai = "s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/human_g1k_v37.fasta.fai"
params.human_dict = "s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/human_g1k_v37.dict"
params.human_fasta = "s3://ucla-rare-diseases/UCLA-UDN/alden_test/rnaseq/references/human_g1k_v37.fasta"

log.info """\
    R N A - S E Q _ W F   P I P E L I N E
    ===================================
    sample_name         : ${params.sample_name}
    library             : ${params.library}
    fastq_bucket        : ${params.fastq_bucket}
    rib_reference_path  : ${params.rib_reference_path}
    outdir              : ${params.outdir}
    """
    .stripIndent(true)

include { download_fastqs; download_rna_ref as download_rrna; download_rna_ref as download_globinrna; download_human_ref } from './modules/download_files.nf'
include { run_fastp } from './modules/fastp.nf'
include { filter_fastq } from './modules/filters.nf'
include { bwa_mem as bwa_mem_rrna; bwa_mem as bwa_mem_globinrna } from './modules/bwa.nf'
include { samtools_view as samtools_view_rrna; samtools_flagstat as samtools_flagstat_rrna; samtools_index; samtools_cram } from './modules/samtools.nf'
include { samtools_view as samtools_view_globinrna; samtools_flagstat as samtools_flagstat_globinrna; samtools_view_sj } from './modules/samtools.nf'
include { check_star_reference; star_alignreads } from './modules/star.nf'
include { bam2sj } from './modules/bam2sj/main.nf'
include { prioritize_splice_junctions } from './modules/prioritize_splice_junctions/main.nf'
include { download_gencode; subread_featurecounts } from './modules/subreads.nf'
include { upload_files } from './modules/upload_outputs.nf'

workflow {
    download_fastqs_ch = download_fastqs(params.sample_name, params.library, params.fastq_bucket)
    download_rrna_ch = download_rrna(params.rib_reference_path, "rrna")
    download_globinrna_ch = download_globinrna(params.rib_reference_path, "globinrna")

    // contamination check
    fastp_ch = run_fastp(download_fastqs_ch)
    filtered_fastq_ch = filter_fastq(fastp_ch)
    rrna_bwa_ch = bwa_mem_rrna(filtered_fastq_ch, download_rrna_ch, "human_rRNA_strict.fasta", "rrna")
    globinrna_bwa_ch = bwa_mem_globinrna(filtered_fastq_ch, download_globinrna_ch, "human_globinRNA.fa", "globinrna")
    rrna_samtools_view_ch = samtools_view_rrna(params.sample_name, rrna_bwa_ch, "rrna")
    globinrna_samtools_view_ch = samtools_view_globinrna(params.sample_name, globinrna_bwa_ch, "globinrna")
    rrna_samtools_flagstat_ch = samtools_flagstat_rrna(params.sample_name, rrna_samtools_view_ch, "rrna")
    globinrna_samtools_flagstat_ch = samtools_flagstat_globinrna(params.sample_name, globinrna_samtools_view_ch, "globinrna")

    // STAR alignment
    star_index_ref_ch = check_star_reference(download_fastqs_ch)
    star_alignreads_ch = star_alignreads(params.sample_name, star_index_ref_ch, fastp_ch)
    samtools_index(star_alignreads_ch)

    // Create SJ tab file
    sam_ch = samtools_view_sj(params.sample_name, star_alignreads_ch) 
    sj_tab_ch = bam2sj(params.sample_name, sam_ch)

    // Create spreadsheet with prioritize splice junctions
    splice_junctions_ch = prioritize_splice_junctions(params.sample_name, sj_tab_ch, params.latest_directory_date, params.latest_udn_id_key, params.sjdblist, params.genemap, params.constraint, params.gencode_cds)

    // Create counts by gene
    gencode_pc_ch = download_gencode(params.gencode_gtf_path)
    feature_counts_ch = subread_featurecounts(params.sample_name, gencode_pc_ch, star_alignreads_ch)

    // Create CRAM files
    download_human_ref_ch = download_human_ref(params.human_fasta, params.human_fai, params.human_dict)
    cram_ch = samtools_cram(params.sample_name, download_human_ref_ch, star_alignreads_ch)

    // Upload selected output files
    upload_files(params.library, rrna_samtools_flagstat_ch, globinrna_samtools_flagstat_ch, star_alignreads_ch, sj_tab_ch, splice_junctions_ch, feature_counts_ch, cram_ch)
}