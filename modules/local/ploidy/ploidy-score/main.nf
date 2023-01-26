process PLOIDY_SCORE {
    tag "ploidy_score"
    label 'process_low'

    container "nicolasvnk/optparse-ubuntu:latest"

    input:
    path ploidy_matrix

    output:
    path "*_ploidy_plots.tar.gz"   , emit: ploidy_plots
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "batch"

    """
    mkdir ploidy_est
    estimatePloidy.R -z -O ./ploidy_est ${ploidy_matrix}

    tar -zcf ./ploidy_est.tar.gz ./ploidy_est
    mv ploidy_est.tar.gz ${prefix}_ploidy_plots.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(Rscript --version 2>&1 | sed -e 's/^R scripting front-end version //; s/ (.*)\$//')
    END_VERSIONS
    """
}
