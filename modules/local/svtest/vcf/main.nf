process SVTEST_VCF {
    tag "$meta.id"
    label 'process_low'

    container "nicolasvnk/svtest:0.1"

    input:
    tuple val(meta), path(vcf), path(tbi), path(baseline_vcf)
    path fasta_fai

    output:
    tuple val(meta), path("*.tsv") , emit: metrics
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = "0.1"

    def arguments   = args.args ?: ''
    def types       = args.types
    def contigs     = fasta_fai

    def baseline = baseline_vcf ? "--baseline-vcf ${baseline_vcf}" : ""

    """
    echo "${meta.id}" > samples.txt

    svtest vcf \\
        ${arguments} \\
        ${baseline} \\
        ${vcf} \\
        ${contigs} \\
        samples.txt \\
        ${types} \\
        ${prefix} \\
        > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svtest: ${VERSION}
        python: \$(python3 --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}
