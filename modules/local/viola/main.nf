process VIOLA {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/nicolasvnk/viola:1.0.2"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_standardized.vcf.gz")  , emit: vcf
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

    variant=\$(cat new_${unzipped_vcf} | awk '/^#/ {next} {print 1;exit}' || echo 0)

    if [ \$variant -eq 1 ]
    then
        viola_standardize.py \\
            new_${unzipped_vcf} \\
            ${meta.caller} \\
            ${prefix}.vcf \\
            ${meta.id}
        bgzip ${prefix}.vcf
        echo "Ran the standardization succesfully"
    else
        echo "${vcf} was empty, so the viola_standardize.py process was skipped."
        bgzip new_${unzipped_vcf}
        cp new_${vcf} ${prefix}.vcf.gz
    fi

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
