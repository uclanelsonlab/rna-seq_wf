
version 1.0

task align_reads {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Task for STAR alignReads"
    }
    input {
        File fastq_r1
        File fastq_r2
        String prefix
        File star_index
        String? star_docker

        # STAR options
        Int? runThreadN
        Int? sjdbOverhang
        Int? limitBAMsortRAM
        Int? outBAMsortingThreadN
    }

    String actual_star_docker=select_first([star_docker, "gvcn/star:2.7.10b"])
    Int actual_runThreadN=select_first([runThreadN, 30])
    Int actual_sjdbOverhang=select_first([sjdbOverhang, 149])
    Int actual_limitBAMsortRAM=select_first([limitBAMsortRAM, 55612616400])
    Int actual_outBAMsortingThreadN=select_first([outBAMsortingThreadN, 16])
    command {
        set -uexo pipefail
        mv ~{star_index} star_index_files.tar.gz
        tar -xvf star_index_files.tar.gz
        mv star_index*/ star_index_files/
        STAR --runMode alignReads \
            --runThreadN ~{actual_runThreadN} \
            --genomeDir $PWD/star_index_files \
            --twopassMode Basic \
            --sjdbOverhang ~{actual_sjdbOverhang} \
            --readFilesIn ~{fastq_r1} ~{fastq_r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix ~{prefix}. \
            --alignSoftClipAtReferenceEnds Yes \
            --quantMode GeneCounts \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMcompression -1 \
            --outSAMunmapped Within \
            --genomeLoad NoSharedMemory \
            --limitBAMsortRAM ~{actual_limitBAMsortRAM} \
            --outBAMsortingThreadN ~{actual_outBAMsortingThreadN} \
            --outSAMattrRGline ID:rg1 SM:~{prefix} PL:Illumina LB:~{prefix}
    }
    runtime {
        docker: actual_star_docker
        dx_instance_type: "mem2_ssd1_v2_x32"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "6H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output { 
        File star_bam = "${prefix}.Aligned.sortedByCoord.out.bam"
        File sjdb_txt = "${prefix}._STARgenome/sjdbInfo.txt"
        File sjdb_tab = "${prefix}._STARgenome/sjdbList.out.tab"
        File read_counts = "${prefix}.ReadsPerGene.out.tab"
        File junctions = "${prefix}.SJ.out.tab"
        File junctions_pass1 = "${prefix}._STARpass1/SJ.out.tab"
        File junctions_pass1_log = "${prefix}._STARpass1/Log.final.out"
        Array[File] logs = ["${prefix}.Log.final.out", "${prefix}.Log.out", "${prefix}.Log.progress.out"]
    }
}

workflow align_reads_wf {
    parameter_meta {
        fastq_r1: {
            description: "Sample R1 FASTQ file",
            extension: ".R1.fastq.gz"
        }
        fastq_r2: {
            description: "Sample R2 FASTQ file",
            extension: ".R2.fastq.gz"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        star_index: {
            description: "TAR zip reference files from Gencode",
            extension: ".tar.gz"
        }
        star_docker: {
            description: "TAR zip docker image from STAR stored on DNAnexus (file-GVjJgfQ02k8bYxZ9z70g869k)",
            extension: ".tar.gz"
        }
    }
    input {
        File fastq_r1
        File fastq_r2
        String prefix
        File star_index
        String? star_docker
    }
    call align_reads {
        input:
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            prefix=prefix,
            star_index=star_index,
            star_docker=star_docker
    }
    output {
        File star_bam = align_reads.star_bam
        File sjdb_txt = align_reads.sjdb_txt
        File sjdb_tab = align_reads.sjdb_tab
        File read_counts = align_reads.read_counts
        File junctions = align_reads.junctions
        File junctions_pass1 = align_reads.junctions_pass1
        File junctions_pass1_log = align_reads.junctions_pass1_log
        Array[File] logs = align_reads.logs
    }
}
