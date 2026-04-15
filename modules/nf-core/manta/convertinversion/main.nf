nextflow.preview.types = true

process MANTA_CONVERTINVERSION {
    tag "${meta.id}"
    label 'process_low'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7b/7b09474b3b6537f95f6fabbd4ed3ae397adad69d195217585e5101c8bdb914aa/data'
        : 'community.wave.seqera.io/library/htslib_manta_samtools_python:0f2533c881652912'}"

    input:
    tuple(meta: Map, vcf: Path)
    tuple(meta2: Map, fasta: Path)

    output:
    vcf = tuple(meta, file("*.vcf.gz"))
    tbi = tuple(meta, file("*.vcf.gz.tbi"))
    versions_manta = tuple("${task.process}", "manta", eval("configManta.py --version"))
    versions_samtools = tuple("${task.process}", "samtools", eval("samtools --version | head -1 | sed -e s'/samtools //'"))

    topic:
    tuple("${task.process}", "manta", eval("configManta.py --version")) >> 'versions'
    tuple("${task.process}", "samtools", eval("samtools --version | head -1 | sed -e s'/samtools //'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convertInversion.py \$(which samtools) ${fasta} ${vcf} | bgzip --threads ${task.cpus} > ${prefix}.vcf.gz
    tabix ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
