nextflow.preview.types = true

process ENSEMBLVEP_DOWNLOAD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3da6e21cbf9803529421d7e136d1ebec5ff71ec50e0d996eda2ce11ec2c19bf9/data'
        : 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf:1e13f65f931a6954'}"

    input:
    tuple(meta: Map, assembly: String, species: String, cache_version: String)

    output:
    cache = tuple(meta, file(prefix))
    versions_ensemblvep = tuple("${task.process}", 'ensemblvep', eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"))
    versions_perlmathcdf = tuple("${task.process}", 'perl-math-cdf', eval("perl -MMath::CDF -e 'print \$Math::CDF::VERSION'"))

    topic:
    tuple("${task.process}", 'ensemblvep', eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'")) >> 'versions'
    tuple("${task.process}", 'perl-math-cdf', eval("perl -MMath::CDF -e 'print \$Math::CDF::VERSION'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'vep_cache'
    """
    vep_install \\
        --CACHEDIR ${prefix} \\
        --SPECIES ${species} \\
        --ASSEMBLY ${assembly} \\
        --CACHE_VERSION ${cache_version} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'vep_cache'
    """
    mkdir ${prefix}
    """
}
