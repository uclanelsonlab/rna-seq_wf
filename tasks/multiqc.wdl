version 1.0

task multiqc_array {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Aggregate results from bioinformatics analyses across many samples into a single report. https://multiqc.info/"
    }
    input {
        Array[File] stats_files
        String prefix
        String? multiqc_image
    }
    String actual_multiqc_image=select_first([multiqc_image, "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"])
    command {
        set -uexo pipefail
        multiqc \
            --file-list ${write_lines(stats_files)} \
            --filename ~{prefix}_multiqc \
            --title ~{prefix} -z
    }
    runtime {
        docker: actual_multiqc_image
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "5H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File multiqc_html = "~{prefix}_multiqc.html"
        File multiqc_data = "~{prefix}_multiqc_data.zip"
    }
}

task multiqc_file {
    meta {
        author: "George Carvalho"
        email: "gcarvalhoneto@mednet.ucla.edu"
        description: "## Aggregate results from bioinformatics analyses across many samples into a single report. https://multiqc.info/"
    }
    input {
        File files_list
        String prefix
        String? multiqc_image
    }
    String actual_multiqc_image=select_first([multiqc_image, "quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0"])
    command {
        set -uexo pipefail
        multiqc \
            --file-list ~{files_list} \
            --filename ~{prefix}_multiqc \
            --title ~{prefix} -z
    }
    runtime {
        docker: actual_multiqc_image
        dx_instance_type: "mem1_ssd1_v2_x16"
        dx_ignore_reuse: true
        dx_restart: object {
            default: 1,
            max: 1,
            errors: object {
                UnresponsiveWorker: 2,
                ExecutionError: 2,
            }
        }
        dx_timeout: "5H30M"
        dx_access: object {
            network: ["*"],
            developer: true
        }
    }
    output {
        File multiqc_html = "~{prefix}_multiqc.html"
        File multiqc_data = "~{prefix}_multiqc_data.zip"
    }
}

workflow multiqc_wf {
    parameter_meta {
        stats_files: {
            description: "All the files created by the QC steps"
        }
        prefix: {
            description: "Sample prefix to be used in output creation"
        }
        multiqc_image: {
            description: "TAR zip docker image from MultiQC stored on DNAnexus (file-GX16qZ802k8xKbkF797xKqP4)"
        }
    }
    input {
        Array[File] stats_files
        String prefix
        String? multiqc_image
    }
    call multiqc_array {
        input:
            stats_files=stats_files,
            prefix=prefix,
            multiqc_image=multiqc_image
    }
    output {
        File multiqc_html = multiqc_array.multiqc_html
        File multiqc_data = multiqc_array.multiqc_data
    }
}