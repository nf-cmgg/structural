nextflow.preview.types = true

process WISECONDORX_PREDICT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0'
        : 'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0'}"

    input:
    tuple(meta: Map, npz: Path)
    tuple(meta2: Map, reference: Path)
    tuple(meta3: Map, blacklist: Path)

    output:
    aberrations_bed = tuple(meta, file("*_aberrations.bed", optional: true))
    bins_bed = tuple(meta, file("*_bins.bed", optional: true))
    segments_bed = tuple(meta, file("*_segments.bed", optional: true))
    chr_statistics = tuple(meta, file("*_statistics.txt", optional: true))
    chr_plots = tuple(meta, file("[!genome_wide]*.png", optional: true))
    genome_plot = tuple(meta, file("genome_wide.png", optional: true))
    versions_wisecondorx = tuple("${task.process}", 'wisecondorx', eval("pip list |& sed -n 's/wisecondorx *//p'"))

    topic:
    tuple("${task.process}", 'wisecondorx', eval("pip list |& sed -n 's/wisecondorx *//p'")) >> 'versions'

    script:
    def args = task.ext.args ?: '--bed --plot'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = blacklist ? "--blacklist ${blacklist}" : ""

    def plots = args.contains("--plot") ? "mv ${prefix}.plots/* ." : ""
    """
    WisecondorX predict \\
        ${npz} \\
        ${reference} \\
        ${prefix} \\
        ${bed} \\
        ${args}

    ${plots}
    """

    stub:
    def args = task.ext.args ?: '--bed --plot'
    def prefix = task.ext.prefix ?: "${meta.id}"

    def bed = args.contains("--bed") ? "touch ${prefix}_aberrations.bed && touch ${prefix}_bins.bed && touch ${prefix}_statistics.txt && touch ${prefix}_segments.bed" : ""
    def plot = args.contains("--plot") ? "touch genome_wide.png && touch chr22.png && touch chr1.png" : ""

    """
    ${bed}
    ${plot}
    """
}
