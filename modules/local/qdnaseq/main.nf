process QDNASEQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-qdnaseq==1.34.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-qdnaseq:1.34.0--r42hdfd78af_0':
        'biocontainers/bioconductor-qdnaseq:1.34.0--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(annotations)

    output:
    tuple val(meta), path("*.bed")              , emit: bed
    tuple val(meta), path("*.cna")              , emit: cna
    tuple val(meta), path("*_segments.txt")     , emit: segments
    tuple val(meta), path("statistics.out")     , emit: statistics
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.34.0"

    """
    # TODO fix the taxonomy ID
    qDNAseq.R \\
        ${bam} \\
        ${prefix} \\
        123456 \\ 
        ${annotations}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qDNAseq: ${VERSION}
    END_VERSIONS
    """

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
