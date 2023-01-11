process ZPASTE {
    tag "zpaste"
    label 'process_low'

    conda 'bioconda::tabix=1.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    path column_files

    output:
    path "*.RD.txt.gz"          , emit: matrix_file
    path "*.RD.txt.gz.tbi"      , emit: matrix_file_index
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "batch"

    """
    # paste unzipped files and compress
    paste ${column_files} | bgzip --threads ${task.cpus} -c > "${prefix}.RD.txt.gz"
    tabix -p bed ${prefix}.RD.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
