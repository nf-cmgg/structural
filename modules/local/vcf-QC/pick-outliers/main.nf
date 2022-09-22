process PICK_OUTLIERS {
    tag "${meta.id}"

    container "nicolasvnk/svtest:0.1"

    input:
    tuple val(meta), path(stat)

    output:
    tuple val(meta), path('*.QC.outlier.low')   , emit: low
    tuple val(meta), path('*.QC.outlier.high')  , emit: high
    path "versions.yml"                         , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    calc_num_svs_pick_outlier.py ${stat} ${prefix}.QC.outlier -z

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}