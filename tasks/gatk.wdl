version 1.0

task gatk_SplitNCigarReads {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK SplitNCigarReads, split reads with N in cigar"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam
        File bai
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.1.3.0"])
    command {
        set -uexo pipefail
        gatk SplitNCigarReads \
            --java-options "-Xmx28G" \
            -R ~{fasta} \
            -I ~{bam} \
            -O ~{prefix}.split_cigar.bam
    }
    runtime {
        docker: actual_gatk_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "8H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File bam_split_cigar = "${prefix}.split_cigar.bam"
        File bai_split_cigar = "${prefix}.split_cigar.bai"
    }
}

task gatk_BaseRecalibrator {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK BaseRecalibrator, it generates recalibration table for Base Quality Score Recalibration (BQSR)"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam
        File bai
        File dbsnp138
        File known_indels
        File indels_1kG
        File af_only_gnomad
        File small_exac_common_3
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.1.3.0"])
    command {
        set -uexo pipefail
        tabix ~{dbsnp138}
        tabix ~{known_indels}
        tabix ~{indels_1kG}
        tabix ~{af_only_gnomad}
        tabix ~{small_exac_common_3}        
        gatk BaseRecalibrator \
            --java-options "-Xmx28G" \
            -R ~{fasta} \
            -I ~{bam} \
            --known-sites ~{dbsnp138} \
            --known-sites ~{known_indels} \
            --known-sites ~{indels_1kG} \
            --known-sites ~{af_only_gnomad} \
            --known-sites ~{small_exac_common_3} \
            -O ~{prefix}.recal_data.table
    }
    runtime {
        docker: actual_gatk_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "8H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File recal_data = "${prefix}.recal_data.table"
    }
}

task gatk_ApplyBQSR {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK ApplyBQSR, apply base quality score recalibration"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam
        File bai
        File recal_data
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.1.3.0"])
    command {
        set -uexo pipefail
        gatk ApplyBQSR \
            --java-options "-Xmx28G" \
            -R ~{fasta} \
            -I ~{bam} \
            --bqsr-recal-file ~{recal_data} \
            -O ~{prefix}.bqsr.bam
    }
    runtime {
        docker: actual_gatk_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "8H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File bqsr_bam = "${prefix}.bqsr.bam"
        File bqsr_bai = "${prefix}.bqsr.bai"
    }
}

task gatk_AnalyzeCovariates {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK AnalyzeCovariates, evaluate and compare base quality score recalibration (BQSR) tables"
    }
    input {
        File recal_data
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.1.3.0"])
    command {
        set -uexo pipefail
        gatk AnalyzeCovariates \
            --java-options "-Xmx28G" \
            -bqsr ~{recal_data} \
            -plots ~{prefix}_AnalyzeCovariates.pdf
    }
    runtime {
        docker: actual_gatk_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "8H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File analyze_pdf = "${prefix}_AnalyzeCovariates.pdf"
    }
}

task gatk_HaplotypeCaller {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK HaplotypeCaller, call germline SNPs and indels via local re-assembly of haplotypes"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam
        File bai
        File dbsnp138
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.1.3.0"])
    command {
        set -uexo pipefail
        tabix ~{dbsnp138}
        gatk HaplotypeCaller \
            --java-options "-Xmx28G" \
            -R ~{fasta} \
            -I ~{bam} \
            -O ~{prefix}_HaplotypeCaller.g.vcf.gz \
            --dbsnp ~{dbsnp138} \
            --add-output-vcf-command-line true \
            -ERC GVCF \
            -G StandardAnnotation \
            -G AS_StandardAnnotation \
            --max-alternate-alleles 3 \
            -contamination 0
    }
    runtime {
        docker: actual_gatk_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "8H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File hc_gvcf = "${prefix}_HaplotypeCaller.g.vcf.gz"
        File hc_gvcf_index = "${prefix}_HaplotypeCaller.g.vcf.gz.tbi"
    }
}

task gatk_GenotypeGVCFs {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK GenotypeGVCFs, perform joint genotyping on one or more samples pre-called with HaplotypeCaller"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File hc_gvcf
        File hc_gvcf_index
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.1.3.0"])
    command {
        set -uexo pipefail
        gatk GenotypeGVCFs \
            --java-options "-Xmx28G" \
            -R ~{fasta} \
            -V ~{hc_gvcf} \
            -O ~{prefix}_HaplotypeCaller.vcf.gz
    }
    runtime {
        docker: actual_gatk_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "8H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File hc_vcf = "${prefix}_HaplotypeCaller.vcf.gz"
        File hc_vcf_index = "${prefix}_HaplotypeCaller.vcf.gz.tbi"
    }
}

task gatk_VariantFiltration {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for GATK VariantFiltration, filter variant calls based on INFO and/or FORMAT annotations"
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File hc_vcf
        File hc_vcf_index
        String prefix
        String? gatk_docker
    }
    String actual_gatk_docker=select_first([gatk_docker, "broadinstitute/gatk:4.1.3.0"])
    command {
        set -uexo pipefail
        gatk VariantFiltration \
            --java-options "-Xmx28G" \
            -R ~{fasta} \
            -V ~{hc_vcf} \
            -O ~{prefix}_VariantFiltration.vcf.gz \
            --filter-name "qc_filter" \
            --filter-expression "AB < 0.2 || MQ0 > 50"
    }
    runtime {
        docker: actual_gatk_docker
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 2,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "8H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File filtered_vcf = "${prefix}_VariantFiltration.vcf.gz"
        File filtered_vcf_index = "${prefix}_VariantFiltration.vcf.gz.tbi"
    }
}

workflow gatk_wf {
    parameter_meta {
        fasta: {
            description: "Unziped reference FASTA file stored by DNAnexus",
            extension: ".fasta"
        }
        fasta_fai: {
            description: "Index for reference FASTA file stored by DNAnexus",
            extension: ".fai"
        }
        fasta_dict: {
            description: "Dictionary for reference FASTA file stored by DNAnexus",
            extension: ".dict"
        }
        bam: {
            description: "BAM aligned by STAR alignReads with vW filter marked duplications",
            extension: ".bam"
        }
        bai: {
            description: "Index for BAM aligned by STAR alignReads with vW filter marked duplications",
            extension: ".bai"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        gatk_docker: {
            description: "GATK docker image (file-GVgyqvQ02k8XGzgY2V8BVVBG)"
        }
    }
    input {
        File fasta
        File fasta_fai
        File fasta_dict
        File bam
        File bai
        File dbsnp138
        File known_indels
        File indels_1kG
        File af_only_gnomad
        File small_exac_common_3
        String prefix
        String? gatk_docker
    }
    call gatk_SplitNCigarReads {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            bam=bam,
            bai=bai,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_BaseRecalibrator {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            bam=gatk_SplitNCigarReads.bam_split_cigar,
            bai=gatk_SplitNCigarReads.bai_split_cigar,
            dbsnp138=dbsnp138,
            known_indels=known_indels,
            indels_1kG=indels_1kG,
            af_only_gnomad=af_only_gnomad,
            small_exac_common_3=small_exac_common_3,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_ApplyBQSR {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            bam=gatk_SplitNCigarReads.bam_split_cigar,
            bai=gatk_SplitNCigarReads.bai_split_cigar,
            recal_data=gatk_BaseRecalibrator.recal_data,
            prefix=prefix,
            gatk_docker=gatk_docker
    }    

    call gatk_AnalyzeCovariates {
        input:
            recal_data=gatk_BaseRecalibrator.recal_data,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_HaplotypeCaller {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            bam=gatk_ApplyBQSR.bqsr_bam,
            bai=gatk_ApplyBQSR.bqsr_bai,
            dbsnp138=dbsnp138,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_GenotypeGVCFs {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            hc_gvcf=gatk_HaplotypeCaller.hc_gvcf,
            hc_gvcf_index=gatk_HaplotypeCaller.hc_gvcf_index,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    call gatk_VariantFiltration {
        input: 
            fasta=fasta,
            fasta_fai=fasta_fai,
            fasta_dict=fasta_dict,
            hc_vcf=gatk_GenotypeGVCFs.hc_vcf,
            hc_vcf_index=gatk_GenotypeGVCFs.hc_vcf_index,
            prefix=prefix,
            gatk_docker=gatk_docker
    }

    output {
        File bqsr_bam = gatk_ApplyBQSR.bqsr_bam
        File bqsr_bai = gatk_ApplyBQSR.bqsr_bai
        File analyze_pdf = gatk_AnalyzeCovariates.analyze_pdf
        File hc_gvcf = gatk_HaplotypeCaller.hc_gvcf
        File hc_gvcf_index = gatk_HaplotypeCaller.hc_gvcf_index
        File hc_vcf = gatk_GenotypeGVCFs.hc_vcf
        File hc_vcf_index = gatk_GenotypeGVCFs.hc_vcf_index
        File filtered_vcf = gatk_VariantFiltration.filtered_vcf
        File filtered_vcf_index = gatk_VariantFiltration.filtered_vcf_index
    }
}
