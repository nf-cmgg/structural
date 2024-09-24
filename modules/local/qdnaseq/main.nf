process QDNASEQ {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/cmgg/qdnaseq:0.0.4"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(annotations)

    output:
    tuple val(meta), path("*.bed")              , emit: bed
    tuple val(meta), path("*.cna")              , emit: cna
    tuple val(meta), path("*_segments.txt")     , emit: segments
    tuple val(meta), path("statistics.out")     , emit: statistics
    path "versions.yml"                         , emit: versions

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
