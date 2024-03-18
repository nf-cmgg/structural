process BCFTOOLS_CONCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':                                                                                                                                                            
        'biocontainers/bcftools:1.18--h8b25389_0' }" 

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    def args2 = task.ext.args2   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"

    def tbi_names = tbis.collect { it.name }

    def tabix_vcfs = vcfs.collect {
        if(tbi_names.contains("${it.name}.tbi" as String)) {
            return ""
        }
        return "    tabix ${it.name}"
    }

    """
    ${tabix_vcfs.join("\n")}
    bcftools concat \\
        $args \\
        --threads $task.cpus \\
        ${vcfs} \\
    | bcftools sort \\
        ${args2} \\
        --output ${prefix}.vcf.gz \\
        --output-type z

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
