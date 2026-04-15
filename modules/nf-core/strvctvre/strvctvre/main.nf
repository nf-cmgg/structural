nextflow.preview.types = true

process STRVCTVRE_STRVCTVRE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/9584eeb6569a511be29d0a07bf80103d59d38715ddb971dddeca0bc72aec41d3/data'
        : 'community.wave.seqera.io/library/liftover_strvctvre:5fec172b808cc48e'}"

    input:
    tuple(meta: Map, sv_file: Path, sv_file_index: Path, assembly: String)
    tuple(meta2: Map, phylop: Path)
    tuple(meta3: Map, data_directory: Path)

    output:
    vcf = tuple(meta, file("*.vcf", optional: true))
    bed = tuple(meta, file("*.bed", optional: true))
    versions_strvctvre = tuple("${task.process}", 'strvctvre', eval("StrVCTVRE.py --help |& sed -n 's/StrVCTVRE: version *//p'"))

    topic:
    tuple("${task.process}", 'strvctvre', eval("StrVCTVRE.py --help |& sed -n 's/StrVCTVRE: version *//p'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = ''
    if (sv_file.name.endsWith('.vcf') || sv_file.name.endsWith('.vcf.gz')) {
        format = 'vcf'
    }
    else if (sv_file.name.endsWith('.bed')) {
        format = 'bed'
    }
    else {
        error("Input structural variants file must be in VCF or BED format")
    }
    if (!['GRCh38', 'GRCh37'].contains(assembly)) {
        error("Assembly must be either 'GRCh37' or 'GRCh38'")
    }
    """
    StrVCTVRE.py \\
        --input ${sv_file} \\
        --format ${format} \\
        --phyloP ${phylop} \\
        --assembly ${assembly} \\
        --liftover liftover_hg19_to_hg38_public.py \\
        --output ${prefix}.${format}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = ''
    if (sv_file.name.endsWith('.vcf') || sv_file.name.endsWith('.vcf.gz')) {
        format = 'vcf'
    }
    else if (sv_file.name.endsWith('.bed')) {
        format = 'bed'
    }
    else {
        error("Input structural variants file must be in VCF or BED format")
    }
    if (!['GRCh38', 'GRCh37'].contains(assembly)) {
        error("Assembly must be either 'GRCh37' or 'GRCh38'")
    }
    """
    touch ${prefix}.${format}
    """
}
