nextflow.preview.types = true

process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/34/34e733a9ae16a27e80fe00f863ea1479c96416017f24a907996126283e7ecd4d/data'
        : 'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b'}"

    input:
    multiqc_files: Set<Path>
    multiqc_config: Path
    extra_multiqc_config: Path
    multiqc_logo: Path
    replace_names: Path
    sample_names: Path

    stage:
    stageAs multiqc_files, "?/*"

    output:
    report = file("*.html")
    data = file("*_data")
    plots = file("*_plots", optional: true)
    versions = tuple("${task.process}", 'multiqc', eval('multiqc --version | sed "s/.* //g"'))

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config ${multiqc_config}" : ''
    def extra_config = extra_multiqc_config ? "--config ${extra_multiqc_config}" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${extra_config} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_report.html
    """
}
