nextflow.enable.dsl = 2

params.sample_name = 'NB-8204-M-muscle'
params.library = "SN_7RNA_S-24-0479_XA044"
params.fastq_bucket = "s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq"
params.rib_reference_path = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference"
params.gencode_gtf_path = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/gencode.v43.primary_assembly.annotation.gtf"
params.outdir = "results"
params.human_fai = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.fa.fai"
params.human_dict = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.dict"
params.human_fasta = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.fa"

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
include { samtools_view as samtools_view_globinrna; samtools_flagstat as samtools_flagstat_globinrna } from './modules/samtools.nf'
include { check_star_reference; star_alignreads } from './modules/star.nf'
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

    // Create counts by gene
    gencode_pc_ch = download_gencode(params.gencode_gtf_path)
    feature_counts_ch = subread_featurecounts(params.sample_name, gencode_pc_ch, star_alignreads_ch)

    // Create CRAM files
    download_human_ref_ch = download_human_ref(params.human_fasta, params.human_fai, params.human_dict)
    cram_ch = samtools_cram(params.sample_name, download_human_ref_ch, star_alignreads_ch)

    // Upload selected output files
    // upload_files(params.library, rrna_samtools_flagstat_ch, globinrna_samtools_flagstat_ch, star_alignreads_ch, feature_counts_ch, cram_ch)
}