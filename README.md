# rna-seq_wf
## RNA-seq pipeline documentation
- Clone the repo
```bash
git clone https://github.com/uclanelsonlab/rna-seq_wf.git
```

- Run the pipeline for your sample, it expects the FASTQ files to be at `s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq` to download
```bash
cd nl-rna-seq_wf/
nextflow run main.nf --sample_name SH1311-P-muscle --library SN_7RNA_S-24-0479_XA044
```

- Check if you have your outputs on S3:
```bash
aws s3 ls s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/ --recursive | grep SH1311-P-muscle
```

- The expected outputs:
    - flagstat_rrna: Ribosomal contamination stats for human_rRNA_strict.fasta (`*.rrna.flagstat.txt`)
    - flagstat_globinrna: Ribosomal contamination stats for human_globinRNA.fa (`*.globinrna.flagstat.txt`)
    - reads_gene: STAR's tab output (`*.ReadsPerGene.out.tab.gz`)
    - reads_gene_log: STAR's log output (`*.ReadsPerGene.log.out`)
    - final_log: STAR's final log (`*.Log.final.out`)
    - sj_tab: STAR's junctions output (`*.SJ.out.tab.gz`)
    - sj_tab_gz: Home made script output of junctions recreated from STAR's output (`*.bam2SJ.out.tab.gz`)
    - all_rare_junctions: Home made script output of all junctions in TSV (`*_rare_junctions_all.tsv`) 
    - rare_junctions: Home made script output of junctions in spreadsheet (`*_rare_junctions_filtered.xlsx`)
    - gene_counts: Subreads featureCounts output (`*.gene_id.exon.ct`)
    - gene_counts_short: Subreads featureCounts output (`*.gene_id.exon.ct.short.txt`)
    - gene_counts_summary: Subreads featureCounts output (`*.gene_id.exon.ct.summary`)
    - rna_cram: CRAM aligned file output from STAR (`*.hg19_rna.normal.cram`)
    - rna_crai: CRAM's index aligned file output from STAR (`*.hg19_rna.normal.cram.crai`)
