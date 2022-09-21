process MAKE_BINCOV_MATRIX_COLUMNS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::tabix=1.11' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(counts_file)
    path binsize
    path bin_locs

    output:
    tuple val(meta), path("*_bincov.RD.txt")   , emit: bincov
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    firstchar=\$(cat ${counts_file} | head -c 1)

    if [ \$firstchar == '@' ]; then
      shift=1  # GATK CollectReadCounts (to convert from 1-based closed intervals)
    else
      shift=0  # bincov sample or matrix
    fi

    binsize=\$(cat ${binsize})

    TMP_BED="${prefix}.tmp.bed"
    printf "#Chr\tStart\tEnd\t%s\n" "${prefix}" > \$TMP_BED
    cat "${counts_file}" \\
      | sed '/^@/d' \\
      | sed '/^CONTIG	START	END	COUNT\$/d' \\
      | sed '/^#/d' \\
      | awk -v x=\$shift -v b=\$binsize \\
        'BEGIN{OFS="\t"}{\$2=\$2-x; if (\$3-\$2==b) print \$0}' \\
      >> "\$TMP_BED"

    cut -f4- "\$TMP_BED" >> "${prefix}_bincov.RD.txt"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
