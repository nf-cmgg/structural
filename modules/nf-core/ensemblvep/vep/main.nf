nextflow.preview.types = true

process ENSEMBLVEP_VEP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3da6e21cbf9803529421d7e136d1ebec5ff71ec50e0d996eda2ce11ec2c19bf9/data'
        : 'community.wave.seqera.io/library/ensembl-vep_perl-math-cdf:1e13f65f931a6954'}"

    input:
    tuple(meta: Map, vcf: Path, custom_extra_files: Path)
    genome: String
    species: String
    cache_version: Integer
    cache: Path
    tuple(meta2: Map, fasta: Path)
    extra_files: Path

    output:
    vcf = tuple(meta, file("${prefix}.vcf.gz", optional: true))
    tbi = tuple(meta, file("${prefix}.vcf.gz.tbi", optional: true))
    tab = tuple(meta, file("${prefix}.tab.gz", optional: true))
    json = tuple(meta, file("${prefix}.json.gz", optional: true))
    report = tuple(meta, "${task.process}", 'ensemblvep', file("*.html", optional: true))
    versions_ensemblvep = tuple("${task.process}", 'ensemblvep', eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'"))
    versions_tabix = tuple("${task.process}", 'tabix', eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"))
    versions_perlmathcdf = tuple("${task.process}", 'perl-math-cdf', eval("perl -MMath::CDF -e 'print \\\$Math::CDF::VERSION'"))

    topic:
    tuple(meta, "${task.process}", 'ensemblvep', file("*.html", optional: true)) >> 'multiqc_files'
    tuple("${task.process}", 'ensemblvep', eval("vep --help | sed -n '/ensembl-vep/s/.*: //p'")) >> 'versions'
    tuple("${task.process}", 'tabix', eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'")) >> 'versions'
    tuple("${task.process}", 'perl-math-cdf', eval("perl -MMath::CDF -e 'print \\\$Math::CDF::VERSION'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json") ? 'json' : args.contains("--tab") ? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    prefix = task.ext.prefix ?: "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def create_index = file_extension == "vcf" ? "tabix ${args2} ${prefix}.${file_extension}.gz" : ""
    """
    vep \\
        -i ${vcf} \\
        -o ${prefix}.${file_extension}.gz \\
        ${args} \\
        ${compress_cmd} \\
        ${reference} \\
        --assembly ${genome} \\
        --species ${species} \\
        --cache \\
        --cache_version ${cache_version} \\
        --dir_cache ${dir_cache} \\
        --fork ${task.cpus}

    ${create_index}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json") ? 'json' : args.contains("--tab") ? 'tab' : 'vcf'
    def create_index = file_extension == "vcf" ? "touch ${prefix}.${file_extension}.gz.tbi" : ""
    """
    echo "" | gzip > ${prefix}.${file_extension}.gz
    ${create_index}
    touch ${prefix}_summary.html
    """
}
