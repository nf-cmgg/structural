process SIMPLE_EVENT_ANNOTATION {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-structuralvariantannotation=1.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-structuralvariantannotation:1.13.0--r42hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-structuralvariantannotation:1.13.0--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: bed
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    simple-event-annotation.R \\
        ${vcf} \\
        ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1 | sed 's/R version //;s/ .*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1 | sed 's/R version //;s/ .*//')
    END_VERSIONS
    """
}
