process ZPASTE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::tabix=1.11' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(counts_file)

    output:
    tuple val(meta), path("") , emit: matrix
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def binsize = task.ext.binsize ?: "NOT_DEFINED"

    """
    # Use named pipes to stream unzipped column files in memory
    mkdir -p column_file_fifos
    FILE_NUM=0
    while read -r COLUMN_FILE; do
      FIFO=$(printf "column_file_fifos/%08d" $FILE_NUM)
      mkfifo "$FIFO"
      bgzip -@$(nproc) -cd "$COLUMN_FILE" > "$FIFO" &
      ((++FILE_NUM))
    done < ~{write_lines(column_files)}

    # paste unzipped files and compress
    paste column_file_fifos/* | bgzip -@$(nproc) -c > "~{matrix_file_name}"
    tabix -p bed "~{matrix_file_name}"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
