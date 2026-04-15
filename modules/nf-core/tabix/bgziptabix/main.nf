nextflow.preview.types = true

process TABIX_BGZIPTABIX {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data'
        : 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'}"

    input:
    tuple(meta: Map, input: Path)

    output:
    gz_index = tuple(meta, file("*.gz"), file("*.{tbi,csi}"))
    versions_tabix = tuple("${task.process}", 'tabix', eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"))
    versions_bgzip = tuple("${task.process}", 'bgzip', eval("bgzip --version | sed '1!d;s/.* //'"))

    topic:
    tuple("${task.process}", 'tabix', eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'")) >> 'versions'
    tuple("${task.process}", 'bgzip', eval("bgzip --version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bgzip --threads ${task.cpus} -c ${args} ${input} > ${prefix}.${input.getExtension()}.gz
    tabix --threads ${task.cpus} ${args2} ${prefix}.${input.getExtension()}.gz

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: ''
    def index = args2.contains("-C ") || args2.contains("--csi") ? "csi" : "tbi"
    """
    echo "" | gzip > ${prefix}.${input.getExtension()}.gz
    touch ${prefix}.${input.getExtension()}.gz.${index}

    """
}
