process SCATTER_BEDS {
    tag "$meta.id"
    label 'process_single'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(bed)
    val(scatter_size)

    output:
    tuple val(meta), path("*.bed")  , emit: scatter
    path('versions.yml')            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    awk -vFS="\t" '{
        if (\$0 ~ /^[^#].*\$/) {
            if (name == "" || size >= ${scatter_size}) {
                name = sprintf("${prefix}_%d.bed", count++)
                size = 0
            }
            size += \$3-\$2+1
            print \$0 > name
        }
    }' ${bed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_1.bed
    touch ${prefix}_2.bed
    touch ${prefix}_3.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
