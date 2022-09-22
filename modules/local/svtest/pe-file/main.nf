process SVTEST_PEFILE {
    tag "$meta.id"
    label 'process_low'

    container "nicolasvnk/svtest:0.1"

    input:
    tuple val(meta), path(pe_file)

    output:
    tuple val(meta), path("*.tsv") , emit: metrics
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = "0.1"
    
    """
    echo "${meta.id}" > samples.txt

    svtest pe-file \\
        ${pe_file} \\
        samples.txt \\
        > ${prefix}.pe-file.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtest: ${VERSION}
        python: \$(python3 --version | sed -e "s/Python v//g")
    END_VERSIONS
    """
}