nextflow.preview.types = true

process MANTA_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f696c93e6209e33ac0d15f1ecfa799bc67329eec07b0569e065ea8b220b53953/data'
        : 'community.wave.seqera.io/library/manta_python:0eb71149179b3920'}"

    input:
    tuple(meta: Map, input: List<Path>, index: List<Path>, target_bed: Path, target_bed_tbi: Path)
    tuple(meta2: Map, fasta: Path)
    tuple(meta3: Map, fai: Path)
    config: Path

    output:
    candidate_small_indels_vcf = tuple(meta, file("*candidate_small_indels.vcf.gz"))
    candidate_small_indels_vcf_tbi = tuple(meta, file("*candidate_small_indels.vcf.gz.tbi"))
    candidate_sv_vcf = tuple(meta, file("*candidate_sv.vcf.gz"))
    candidate_sv_vcf_tbi = tuple(meta, file("*candidate_sv.vcf.gz.tbi"))
    diploid_sv_vcf = tuple(meta, file("*diploid_sv.vcf.gz"))
    diploid_sv_vcf_tbi = tuple(meta, file("*diploid_sv.vcf.gz.tbi"))
    versions_manta = tuple("${task.process}", "manta", eval("configManta.py --version"))

    topic:
    tuple("${task.process}", "manta", eval("configManta.py --version")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_files = input.collect { bam -> "--bam ${bam}" }.join(' ')
    def options_manta = target_bed ? "--callRegions ${target_bed}" : ""
    def config_option = config ? "--config ${config}" : ""
    """
    configManta.py \\
        ${input_files} \\
        ${config_option} \\
        --reference ${fasta} \\
        --runDir manta \\
        ${options_manta} \\
        ${args}

    python manta/runWorkflow.py -m local -j ${task.cpus}

    mv manta/results/variants/candidateSmallIndels.vcf.gz \\
        ${prefix}.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
        ${prefix}.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \\
        ${prefix}.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \\
        ${prefix}.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz \\
        ${prefix}.diploid_sv.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi \\
        ${prefix}.diploid_sv.vcf.gz.tbi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.candidate_small_indels.vcf.gz
    touch ${prefix}.candidate_small_indels.vcf.gz.tbi
    echo "" | gzip > ${prefix}.candidate_sv.vcf.gz
    touch ${prefix}.candidate_sv.vcf.gz.tbi
    echo "" | gzip > ${prefix}.diploid_sv.vcf.gz
    touch ${prefix}.diploid_sv.vcf.gz.tbi
    """
}
