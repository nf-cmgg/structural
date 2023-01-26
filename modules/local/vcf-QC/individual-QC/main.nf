process INDVIDUAL_QC {
    tag "${meta.id}"

    container "nicolasvnk/svtest:0.1"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.QC.stat')   , emit: stat
    path "versions.yml"                  , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    calcu_num_SVs.by_type_chromo.py ${vcf} ${prefix}.QC.stat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
