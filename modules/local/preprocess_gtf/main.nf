process PREPROCESS_GTF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8a/8ad257d53c2a2b8810d2b12d4d8e3ea438bc8c4a6be7c39b0354cd7bb8d5c260/data':
        'community.wave.seqera.io/library/python:3.13.5--18032a8dc5d4b91e' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.sanitized.gtf"), emit: gtf
    tuple val("${task.process}"), val('python'), eval("python --version |& sed '1!d; s/^Python //g'"), emit: versions_python, topic: versions

    script:
    def prefix  = task.ext.prefix ?: "${gtf.baseName}"

    """
    preprocess_gtf.py $gtf ${prefix}.sanitized.gtf
    """

    stub:
    def prefix = task.ext.prefix ?: "${gtf.baseName}"

    """
    touch ${prefix}.sanitized.gtf
    """
}
