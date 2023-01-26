process CALCMEDCOV {
    tag "calcmedcov"
    label 'process_low'

    container "nicolasvnk/optparse-ubuntu:latest"

    input:
    path bincov_matrix

    output:
    path "*_medianCov.transposed.bed"   , emit: median_cov_file
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "batch"

    """
    zcat ${bincov_matrix} > ${prefix}_fixed.bed
    medianCoverage.R ${prefix}_fixed.bed -H ${prefix}_medianCov.bed
    Rscript -e "x <- read.table(\"${prefix}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"${prefix}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(Rscript --version 2>&1 | sed -e 's/^R scripting front-end version //; s/ (.*)\$//')
    END_VERSIONS
    """
}
