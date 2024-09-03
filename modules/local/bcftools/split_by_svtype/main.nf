process BCFTOOLS_SPLIT_BY_SVTYPE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: split_vcfs
    tuple val(meta), path("*.txt")                    , emit: header
    path "versions.yml"                               , emit: versions

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"

    def types = ["del", "ins", "inv", "tra", "bnd", "dup", "other"]
    def check_variants = types.collect { type ->
        def type_vcf = "${prefix}.${type}.${extension}"
        return """
        if [[ -z \$(bcftools view -H ${type_vcf}) ]]; then
            rm ${type_vcf}
        fi

        """
    }

    """
    bcftools view ${args} --threads ${task.cpus} --include 'INFO/SVTYPE="DEL"' --output ${prefix}.del.${extension} ${vcf}
    bcftools view ${args} --threads ${task.cpus} --include 'INFO/SVTYPE="INS"' --output ${prefix}.ins.${extension} ${vcf}
    bcftools view ${args} --threads ${task.cpus} --include 'INFO/SVTYPE="INV"' --output ${prefix}.inv.${extension} ${vcf}
    bcftools view ${args} --threads ${task.cpus} --include 'INFO/SVTYPE="TRA"' --output ${prefix}.tra.${extension} ${vcf}
    bcftools view ${args} --threads ${task.cpus} --include 'INFO/SVTYPE="BND"' --output ${prefix}.bnd.${extension} ${vcf}
    bcftools view ${args} --threads ${task.cpus} --include 'INFO/SVTYPE="DUP"' --output ${prefix}.dup.${extension} ${vcf}
    bcftools view ${args} --threads ${task.cpus} --exclude 'INFO/SVTYPE="DEL" || INFO/SVTYPE="INS" || INFO/SVTYPE="INV" || INFO/SVTYPE="TRA" || INFO/SVTYPE="BND" || INFO/SVTYPE="DUP"' --output ${prefix}.other.${extension} ${vcf}

    bcftools view -h ${vcf} | grep -E '(##INFO|##FORMAT|##ALT|##FILTER)' > ${prefix}.header.txt

    ${check_variants.join("\n")}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"
    """
    echo "" | gzip > ${prefix}.del.${extension}
    echo "" | gzip > ${prefix}.ins.${extension}
    echo "" | gzip > ${prefix}.inv.${extension}
    echo "" | gzip > ${prefix}.tra.${extension}
    echo "" | gzip > ${prefix}.bnd.${extension}
    echo "" | gzip > ${prefix}.dup.${extension}
    echo "" | gzip > ${prefix}.other.${extension}

    touch ${prefix}.header.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
