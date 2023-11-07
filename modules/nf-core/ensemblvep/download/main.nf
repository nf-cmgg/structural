process ENSEMBLVEP_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    conda "${ DynamicContainers.findLatest('ensembl-vep', params.vep_version, 'conda') }"
    container "${ DynamicContainers.findLatest('ensembl-vep', params.vep_version, workflow.containerEngine) }"

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)

    output:
    tuple val(meta), path("vep_cache"), emit: cache
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vep_install \\
        --CACHEDIR vep_cache \\
        --SPECIES $species \\
        --ASSEMBLY $assembly \\
        --CACHE_VERSION $cache_version \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir vep_cache

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
