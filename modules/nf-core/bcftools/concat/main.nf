process BCFTOOLS_CONCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"

    def tbi_names = tbis.collect { it.name }

    def tabix_vcfs = vcfs.collect {
        if(tbi_names.contains("${it.name}.tbi" as String)) {
            return ""
        }
        return "tabix ${it.name}"
    }

    """
    ${tabix_vcfs.join("\n")}
    bcftools concat \\
        --output ${prefix}.vcf.gz \\
        $args \\
        --threads $task.cpus \\
        ${vcfs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
