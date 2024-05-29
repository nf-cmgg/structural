process ANNOTSV_ANNOTSV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/annotsv:3.3.6--py311hdfd78af_0' :
        'biocontainers/annotsv:3.3.6--py311hdfd78af_0' }"

    input:
    tuple val(meta), path(sv_vcf), path(sv_vcf_index), path(candidate_small_variants)
    tuple val(meta2), path(annotations)
    tuple val(meta3), path(candidate_genes)
    tuple val(meta4), path(false_positive_snv)
    tuple val(meta5), path(gene_transcripts)

    output:
    tuple val(meta), path("*.tsv")              , emit: tsv
    tuple val(meta), path("*.unannotated.tsv")  , emit: unannotated_tsv, optional: true
    tuple val(meta), path("*.vcf")              , emit: vcf
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def cand_genes = candidate_genes ? "-candidateGenesFile ${candidate_genes}" : ""
    def small_variants = candidate_small_variants ? "-candidateSnvIndelFiles ${candidate_small_variants}" : ""
    def fp_snv = false_positive_snv ? "-snvIndelFiles ${false_positive_snv}" : ""
    def transcripts = gene_transcripts ? "-txFile ${gene_transcripts}" : ""

    """
    AnnotSV \\
        -annotationsDir ${annotations} \\
        ${cand_genes} \\
        ${small_variants} \\
        ${fp_snv} \\
        ${transcripts} \\
        -outputFile ${prefix}.raw.tsv \\
        -SVinputFile ${sv_vcf} \\
        ${args}

    mv *_AnnotSV/* .
    awk 'BEGIN { FS=OFS="\t" } { if (NR > 1 && NF >= 8) \$1 = \$1 "_" NR; print }' ${prefix}.raw.tsv > ${prefix}.tsv

    variantconvert convert -i ${prefix}.tsv -o ${prefix}.vcf -fi annotsv -fo vcf -c /usr/local/share/python3/variantconvert/configs/GRCh38/annotsv3_from_vcf.json
    sed -i 's/contig=<ID=MT/contig=<ID=M/' ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help 2>&1 | head -n1 | sed 's/^AnnotSV //'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    touch ${prefix}.unannotated.tsv
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help 2>&1 | head -n1 | sed 's/^AnnotSV //'))
    END_VERSIONS
    """
}
