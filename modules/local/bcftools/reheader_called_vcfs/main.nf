process REHEADER_CALLED_VCFS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1':
        'biocontainers/bcftools:1.16--hfe4b78e_1' }"

    input:
    tuple val(meta), path(vcf)
    path header
    path fai

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def update_sequences = fai ? "-f $fai" : ""

    """
    bcftools view -h ${vcf} | grep -E \\#\\#fileformat=VCF\\|\\#\\#fileDate=\\|\\#\\#reference=\\|\\#\\#smoove_counts_stats= \\
        > new_header.txt

    cat ${header} >> new_header.txt

    bcftools view -h ${vcf} | grep -E \\#CHROM >> new_header.txt

    bcftools \\
        reheader \\
        ${update_sequences} \\
        -h new_header.txt \\
        ${args} \\
        --threads ${task.cpus} \\
        -o ${prefix}.vcf.gz \\
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
