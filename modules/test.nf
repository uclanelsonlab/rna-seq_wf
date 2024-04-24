process echo_files {
    tag "Echo files $meta"

    input:
    tuple val(meta), path(reads)

    output:
    stdout

    script:

    """
    echo ${meta} ${reads[0]} ${reads[1]}
    """
}