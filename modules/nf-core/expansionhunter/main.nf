nextflow.preview.types = true

process EXPANSIONHUNTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/14e1d96665f934a98e569fc5a6fa237f98d3753eee2b6f60d0aea8ff9d44f406/data'
        : 'community.wave.seqera.io/library/expansionhunter:5.0.0--389ada7e191a4fba'}"

    input:
    tuple(meta: Map, bam: Path, bai: Path)
    tuple(meta2: Map, fasta: Path)
    tuple(meta3: Map, fasta_fai: Path)
    tuple(meta4: Map, variant_catalog: Path)

    output:
    vcf = tuple(meta, file("*.vcf.gz"))
    tbi = tuple(meta, file("*.vcf.gz.tbi"))
    json = tuple(meta, file("*.json.gz"))
    bam = tuple(meta, file("*_realigned.bam"))
    versions_expansionhunter = tuple("${task.process}", 'expansionhunter', eval("ExpansionHunter --version | head -1 | sed -n 's/^.*ExpansionHunter v//; s/]//p'"))
    versions_bgzip = tuple("${task.process}", 'bgzip', eval("bgzip --version | sed '1!d;s/.* //'"))

    topic:
    tuple("${task.process}", 'expansionhunter', eval("ExpansionHunter --version | head -1 | sed -n 's/^.*ExpansionHunter v//; s/]//p'")) >> 'versions'
    tuple("${task.process}", 'bgzip', eval("bgzip --version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ExpansionHunter \\
        ${args} \\
        --reads ${bam} \\
        --output-prefix ${prefix} \\
        --reference ${fasta} \\
        --variant-catalog ${variant_catalog}

    bgzip --threads ${task.cpus} ${args2} ${prefix}.vcf
    tabix ${prefix}.vcf.gz
    bgzip --threads ${task.cpus} ${args2} ${prefix}.json

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.json.gz
    touch ${prefix}_realigned.bam

    """
}
