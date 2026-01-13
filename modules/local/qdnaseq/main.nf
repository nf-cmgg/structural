process QDNASEQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/77b057272e6af69070dc1ee73d0a39d144e6641c0f4e625673de979b21b7bfd0/data':
        'community.wave.seqera.io/library/bioconductor-qdnaseq_r-base_r-lsr:0304e1e0cbed3eab' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(annotations)

    output:
    tuple val(meta), path("*.bed")              , emit: bed
    tuple val(meta), path("*.cna")              , emit: cna
    tuple val(meta), path("*_segments.txt")     , emit: segments
    tuple val(meta), path("statistics.out")     , emit: statistics
    path "versions.yml"                         , emit: versions, topic: versions

    script:
    template "qDNAseq.R"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.34.0"

    """
    touch ${prefix}.bed
    touch ${prefix}.cna
    touch ${prefix}_segments.txt
    touch statistics.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qDNAseq: ${VERSION}
    END_VERSIONS
    """
}
