process BCFTOOLS_CONSENSUS_REHEADER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(merged_vcf), path(single_vcfs), path(single_tbis)
    tuple val(meta2), path(fai)
    val(additional_headers)

    output:
    tuple val(meta), path("*.${extension}"), emit: vcf
    path "versions.yml"                    , emit: versions

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def fai_argument      = fai ? "--fai $fai" : ""
    def add_additional = additional_headers ?
    """
    cat <<-EOF >> ${prefix}.temp.txt
${additional_headers.join("\t\t\n")}
    EOF
    """ : ""

    def args2 = task.ext.args2 ?: '--output-type z'
    extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                    args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                    args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                    args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                    "vcf"
    """
    touch ${prefix}.temp.txt
    for FILE in ${merged_vcf} ${single_vcfs};
    do
        bcftools view --threads ${task.cpus} -h \$FILE | grep -vE '^(#CHROM|##fileformat|##filedate|##contig)' >> ${prefix}.temp.txt
    done

    # Add fileformat and file date from the merged VCF
    bcftools view --threads ${task.cpus} -h ${merged_vcf} \\
        | grep -E '(##filedate|##fileformat)' \\
        > ${prefix}.header.vcf

    # Add aditional header lines
    ${add_additional}

    # Remove duplicate header fields, sort and add them to the header
    awk -F, '!seen[substr(\$0, index(\$0, "##") + 1, index(\$0, ",") - index(\$0, "##") - 1)]++' ${prefix}.temp.txt \\
        | sort \\
        >> ${prefix}.header.vcf

    # Add the CHROM and contig lines to the header
    bcftools view --threads ${task.cpus} -h ${merged_vcf} \\
        | grep -E '(#CHROM|##contig)' \\
        >> ${prefix}.header.vcf

    bcftools \\
        reheader \\
        $fai_argument \\
        --header ${prefix}.header.vcf \\
        $args \\
        --threads $task.cpus \\
        $merged_vcf \\
        | bcftools view \\
        $args2 \\
        --output ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    def args2 = task.ext.args2 ?: '--output-type z'
    def prefix = task.ext.prefix ?: "${meta.id}"

    extension = args2.contains("--output-type b") || args2.contains("-Ob") ? "bcf.gz" :
                    args2.contains("--output-type u") || args2.contains("-Ou") ? "bcf" :
                    args2.contains("--output-type z") || args2.contains("-Oz") ? "vcf.gz" :
                    args2.contains("--output-type v") || args2.contains("-Ov") ? "vcf" :
                    "vcf"
    """
    echo "" | gzip > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
