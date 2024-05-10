# rna-seq_wf
## RNA-seq pipeline documentation
- Clone the repo
```bash
git clone https://github.com/uclanelsonlab/rna-seq_wf.git
```

- Run the pipeline for your sample, it expects the FASTQ files to be at `s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq` to download
```bash
cd nl-rna-seq_wf/
nextflow run main.nf --sample_name SH1311-P-muscle --library SN_7RNA_S-24-0479_XA044 -with-report SH1311-P-muscle.html
```

- Check if you have your outputs on S3:
```bash
aws s3 ls s3://ucla-rare-diseases/UCLA-UDN/rnaseq/output/ --recursive | grep SH1311-P-muscle
```

- The expected outputs:
flagstat_rrna: 
flagstat_globinrna: 
reads_gene: 
reads_gene_log: 
final_log: 
sj_tab: 
sj_tab_gz: 
all_rare_junctions: 
rare_junctions: 
gene_counts:
gene_counts_short: 
gene_counts_summary: 
rna_cram: 
rna_crai: 
