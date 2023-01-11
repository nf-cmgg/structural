process MAKE_BINCOV_MATRIX_COLUMNS {
    tag "$meta.id"
    label 'process_low'

    conda 'bioconda::tabix=1.11'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(count_file)
    path bin_locs

    output:
    tuple val(meta), path("*_bincov.RD.txt.gz")    , emit: bincov
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def binsize = meta.binsize

    """
    firstchar=\$(head -c 1 ${count_file})

    if [ \$firstchar == '@' ]; then
      shift=1  # GATK CollectReadCounts (to convert from 1-based closed intervals)
    else
      shift=0  # bincov sample or matrix
    fi

    # Turn the counts file into a 0-based interval file (BED format)
    TMP_BED="${prefix}.tmp.bed"
    printf "#Chr\\tStart\\tEnd\\t${prefix}\\n" > \$TMP_BED
    cat "${count_file}" \\
      | sed '/^@/d' \\
      | sed '/^CONTIG	START	END	COUNT\$/d' \\
      | sed '/^#/d' \\
      | awk -v x=\$shift \\
        'BEGIN{OFS="\\t"}{\$2=\$2-x; if (\$3-\$2==${binsize}) print \$0}' \\
      >> "\$TMP_BED"

    if ! cut -f1-3 "\$TMP_BED" | cmp <(bgzip -cd ${bin_locs}); then
      echo "${count_file} has different intervals than ${bin_locs}"
      exit 1
    fi

    cut -f4- "\$TMP_BED" | bgzip -c >> ${prefix}.RD.txt.gz



    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
END_VERSIONS
    """
}
