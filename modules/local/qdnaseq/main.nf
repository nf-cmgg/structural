nextflow.preview.types = true

process QDNASEQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/77b057272e6af69070dc1ee73d0a39d144e6641c0f4e625673de979b21b7bfd0/data'
        : 'community.wave.seqera.io/library/bioconductor-qdnaseq_r-base_r-lsr:0304e1e0cbed3eab'}"

    input:
    tuple(meta: Map, bam: Path, bai: Path)
    tuple(meta2: Map, annotations: Path)

    output:
    bed = tuple(meta, file("*.bed"))
    cna = tuple(meta, file("*.cna"))
    segments = tuple(meta, file("*_segments.txt"))
    statistics = tuple(meta, file("statistics.out"))

    topic:
    file("versions.yml") >> 'versions'

    script:
    template("qDNAseq.R")

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
