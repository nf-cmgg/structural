process REVERSE_BED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bed)
    path fasta_fai

    output:
    tuple val(meta), path("*_reversed.bed") , emit: bed
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk 'BEGIN {FS="\\t"}; {print \$1 FS "0" FS \$2}' ${fasta_fai} > ${fasta_fai}.bed

    bedtools \\
        subtract \\
        -a ${fasta_fai}.bed \\
        -b ${bed} \\
        ${args} \\
        > ${prefix}_reversed.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
