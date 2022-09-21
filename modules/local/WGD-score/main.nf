process WGD_SCORE {
    tag "wgd_score"
    label 'process_low'

    container "nicolasvnk/optparse-ubuntu:latest"

    input:
    tuple path(wgd_matrix), path(wgd_scoring_mask) 

    output:
    path "*_WGD_scores.txt.gz"          , emit: wgd_scores
    path "*_WGD_score_distributions.pdf", emit: wgd_scores
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "batch"

    """
    scoreDosageBiases.R -z -O . ${wgd_matrix} ${wgd_scoring_mask}
    mv WGD_scores.txt.gz ${prefix}_WGD_scores.txt.gz
    mv WGD_score_distributions.pdf ${prefix}_WGD_score_distributions.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(Rscript --version 2>&1 | sed -e 's/^R scripting front-end version //; s/ (.*)\$//')
    END_VERSIONS
    """
}
