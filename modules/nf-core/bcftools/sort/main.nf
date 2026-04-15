nextflow.preview.types = true

process BCFTOOLS_SORT {
    tag "${input.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    input: BcftoolsSortInput

    output:
    input + record(
        vcf: file("*.{vcf,vcf.gz,bcf,bcf.gz}"),
        tbi: file("*.tbi", optional: true),
        csi: file("*.csi", optional: true)
    )

    topic:
    tuple("${task.process}", 'bcftools', eval("bcftools --version | sed '1!d; s/^.*bcftools //'")) >> 'versions'

    script:
    def args = task.ext.args ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${input.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def max_memory = task.memory ? "--max-mem ${task.memory.toUnit('MB') * 0.9}M" : ""
    """
    bcftools \\
        sort \\
        --output ${prefix}.${extension} \\
        --temp-dir . \\
        ${max_memory} \\
        ${args} \\
        ${input.vcf}
    """

    stub:
    def args = task.ext.args ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${input.id}"

    def extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def index = args.contains("--write-index=tbi") || args.contains("-W=tbi")
        ? "tbi"
        : args.contains("--write-index=csi") || args.contains("-W=csi")
            ? "csi"
            : args.contains("--write-index") || args.contains("-W")
                ? "csi"
                : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index ==~ "csi|tbi" ? "touch ${prefix}.${extension}.${index}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}

record BcftoolsSortInput {
    id: String
    vcf: Path
}