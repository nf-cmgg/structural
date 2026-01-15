process FIX_CALLERS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/41b4cdbfbbd80bb6f5eea8e93b84b6f2c9c4fb9f64293eb36c7e2e2ae319a062/data':
        'community.wave.seqera.io/library/htslib_python:94ea753cc6037594' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz")  , emit: vcf
    path "versions.yml"                , emit: versions

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
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed s'/Python //')
        bgzip: \$(bgzip --version | head -1 | sed s'/bgzip (htslib) //')
    END_VERSIONS
    """
}
