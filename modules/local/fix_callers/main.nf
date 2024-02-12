process FIX_CALLERS {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/cmgg/python-tabix:0.0.1"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz")  , emit: vcf
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (["${prefix}.vcf.gz", "${prefix}.vcf"].contains("${vcf}")) {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    """
    fix_callers.py $vcf ${prefix}.vcf
    bgzip ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed s'/Python //')
        bgzip: \$(bgzip --version | head -1 | sed s'/bgzip (htslib) //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (["${prefix}.vcf.gz", "${prefix}.vcf"].contains("${vcf}")) {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed s'/Python //')
        bgzip: \$(bgzip --version | head -1 | sed s'/bgzip (htslib) //')
    END_VERSIONS
    """
}
