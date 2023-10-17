process REHEADER_CALLED_VCFS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf)
    path header
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.${extension}")   , emit: vcf
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def update_sequences = fai ? "-f $fai" : ""

    extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                "vcf"

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
        ${vcf} \\
        | bcftools convert ${args2} --output ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: ''

    extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                "vcf"

    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
