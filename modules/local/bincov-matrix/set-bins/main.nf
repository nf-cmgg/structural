process SET_BINS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::tabix=1.11' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(counts_file)

    output:
    path "locs.bed.gz"      , emit: bin_locs
    path "binsize.txt"      , emit: binsize
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def binsize = task.ext.binsize ?: "NOT_DEFINED"

    """
    # make the CollectReadCounts output consistent with the old bincov code
    # determine what format this is

    firstchar=\$(head -c 1 ${counts_file})

    if [ \$firstchar == '@' ]; then
      shift=1  # GATK CollectReadCounts (to convert from 1-based closed intervals)
    else
      shift=0  # bincov sample or matrix
    fi

    # kill the dictionary | kill the header | adjust to bed format: 0-based half-open intervals
    cat ${counts_file} \\
      | sed '/^@/d' \\
      | sed '/^CONTIG	START	END	COUNT\$/d' \\
      | sed '/^#/d' \\
      | awk -v x="\${shift}" 'BEGIN{OFS="\t"}{\$2=\$2-x; print \$1,\$2,\$3}' > tmp_locs

    # determine bin size, and drop all bins not exactly equal to this size
    if [ ${binsize} == 'NOT_DEFINED' ]; then
      # use the most common bin size from the bins
      binsize=\$(
        sed -n '1,1000p' tmp_locs | awk '{ print \$3-\$2 }' \\
        | sort | uniq -c | sort -nrk1,1 \\
        | sed -n '1p' | awk '{ print \$2 }'
      )
    else
      # use the provided bin size
      binsize=${binsize}
    fi

    # store binsize
    echo \$binsize > binsize.txt

    # write final bed file with header, and compress it
    awk -v FS="\t" -v b=\$binsize 'BEGIN{ print "#Chr\tStart\tEnd" } { if (\$3-\$2==b) print \$0 }' tmp_locs \\
        | bgzip -c \\
        > "locs.bed.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
