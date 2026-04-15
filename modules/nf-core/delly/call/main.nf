nextflow.preview.types = true

process DELLY_CALL {
    tag "${input.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/delly:1.3.3--h4d20210_0'
        : 'biocontainers/delly:1.3.3--h4d20210_0'}"

    input:
    input: DellyCallInput

    output:
    input + record(bcf: file("*.{bcf,vcf.gz}"), csi: file("*.{csi,tbi}"))

    topic:
    tuple("${task.process}", 'delly', eval("delly --version |& sed -n '1s/Delly version: *v//p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${input.id}"
    def suffix = task.ext.suffix ?: "bcf"

    def exclude = input.exclude_bed ? "--exclude ${input.exclude_bed}" : ""

    def bcf_output = suffix == "bcf" ? "--outfile ${prefix}.bcf" : ""
    def vcf_output = suffix == "vcf" ? "| bgzip ${args2} --threads ${task.cpus} --stdout > ${prefix}.vcf.gz && tabix ${prefix}.vcf.gz" : ""

    def genotype = input.vcf ? "--vcffile ${input.vcf}" : ""

    """
    delly \\
        call \\
        ${args} \\
        ${bcf_output} \\
        --genome ${input.fasta} \\
        ${genotype} \\
        ${exclude} \\
        ${input.input} \\
        ${vcf_output}
    """

    stub:
    def prefix = task.ext.prefix ?: "${input.id}"
    def suffix = task.ext.suffix ?: "bcf"

    def bcf_output = suffix == "bcf" ? "touch ${prefix}.bcf && touch ${prefix}.bcf.csi" : ""
    def vcf_output = suffix == "vcf" ? "echo '' | gzip > ${prefix}.vcf.gz && touch ${prefix}.vcf.gz.tbi" : ""

    """
    ${bcf_output}
    ${vcf_output}
    """
}

record DellyCallInput {
    id: String
    input: Path
    input_index: Path
    fasta: Path
    fai: Path
    vcf: Path?
    vcf_index: Path?
    exclude_bed: Path?
}
