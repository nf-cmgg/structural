process VIOLA {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/nicolasvnk/viola:1.0.2"

    input:
    tuple val(meta), path(vcf)
    val(caller)

    output:
    tuple val(meta), path("*.vcf.gz")  , emit: vcf
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0.2"

    if ("${vcf}" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }

    def unzipped_vcf = vcf.name.replace(".gz","")

    """
    cp ${vcf} new_${vcf}
    bgzip -d new_${vcf}

    viola_standardize.py \\
        new_${unzipped_vcf} \\
        ${caller} \\
        ${prefix}.vcf \\
        ${meta.id}
    bgzip ${prefix}.vcf
    echo "Ran the standardization succesfully"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viola-sv: ${VERSION}
        python: \$(python3 --version | sed -e "s/Python //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0.2"

    if ("${vcf}" == "${prefix}.vcf.gz") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }


    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        viola-sv: ${VERSION}
        python: \$(python3 --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}
