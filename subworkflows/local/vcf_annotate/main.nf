include { VCF_ANNOTATE_VCFANNO  } from '../../../subworkflows/local/vcf_annotate_vcfanno/main'

include { ENSEMBLVEP_VEP        } from '../../../modules/nf-core/ensemblvep/vep/main'
include { GATK4_SVANNOTATE      } from '../../../modules/nf-core/gatk4/svannotate/main'
include { STRVCTVRE_STRVCTVRE   } from '../../../modules/nf-core/strvctvre/strvctvre/main'
include { TABIX_BGZIPTABIX      } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow VCF_ANNOTATE {
    take:
    vcfs                                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
    fasta                                // channel: [mandatory] [ val(meta), path(fasta) ] => The fasta reference file
    fai                                  // channel: [optional]  [ val(meta), path(fai) ] => The fasta index file
    dict                                 // channel: [optional]  [ val(meta), path(dict) ] => The fasta dict file
    gtf                                  // channel: [optional]  [ val(meta), path(gtf) ] => The preprocessed GTF file for SVAnnotate
    vep_cache                            // channel: [optional]  [ path(cache) ] => The path to the local VEP cache
    vep_extra_files                      // channel: [optional]  [ path(file1, file2, file3...) ] => The VEP extra files
    vcfanno_lua                          // channel: [optional]  [ path(lua) ] => The lua script to influence VCFanno
    vcfanno_toml                         // channel: [optional]  [ path(toml) ] => A vcfanno TOML config
    strvctvre_phylop                     // channel: [optional]  [ val(meta), path(phylop) ] => The phylop bigwig file for StrVCTVRE
    strvctvre_data                       // channel: [optional]  [ val(meta), path(data_dir) ] => The data directory for StrVCTVRE
    genome                               // string:  [mandatory] => The genome used by the variant callers
    species                              // string:  [mandatory] => The species used by VEP
    vep_cache_version                    // integer: [mandatory] => The VEP cache version to use
    vcfanno_resources                    // list:    [optional]  [ path(file1, file2, file3...) ] => The extra VCFanno files
    vcfanno_default_tomls                // list:    [mandatory] => A list of default vcfanno configs to be concatenated with the input TOML
    tools                                // list:    [optional]  => The tools used for annotation (default: all tools)

    main:
    def ch_versions = Channel.empty()
    def ch_reports = Channel.empty()

    def ch_vep = vcfs
    if( tools.contains("vep") || tools.contains("all") ) {
        ENSEMBLVEP_VEP(
            vcfs.map { meta, vcf, _tbi -> [ meta, vcf, []]},
            genome,
            species,
            vep_cache_version,
            vep_cache,
            fasta,
            vep_extra_files
        )
        ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())
        ch_reports  = ch_reports.mix(ENSEMBLVEP_VEP.out.report)
        ch_vep = ENSEMBLVEP_VEP.out.vcf.join(ENSEMBLVEP_VEP.out.tbi, failOnDuplicate:true, failOnMismatch:true)
    }

    def ch_vcfanno = ch_vep
    if( tools.contains("vcfanno") || tools.contains("all") ) {
        VCF_ANNOTATE_VCFANNO(
            ch_vep.map { meta, vcf, tbi -> [ meta, vcf, tbi, []]},
            vcfanno_lua,
            vcfanno_resources,
            vcfanno_toml,
            vcfanno_default_tomls
        )
        ch_versions = ch_versions.mix(VCF_ANNOTATE_VCFANNO.out.versions)
        ch_vcfanno = VCF_ANNOTATE_VCFANNO.out.vcfs
    }

    def ch_svannotate = ch_vcfanno
    if( tools.contains("svannotate") || tools.contains("all") ) {
        def ch_svannotate_input = ch_vcfanno
            .map { meta, vcf, tbi ->
                [ meta, vcf, tbi, [], [] ] // TODO add BED files
            }

        GATK4_SVANNOTATE(
            ch_svannotate_input,
            fasta,
            fai,
            dict,
            gtf
        )
        ch_versions = ch_versions.mix(GATK4_SVANNOTATE.out.versions.first())
        ch_svannotate = GATK4_SVANNOTATE.out.vcf.join(GATK4_SVANNOTATE.out.tbi, failOnMismatch:true, failOnDuplicate:true)
    }

    def ch_strvctvre = ch_svannotate
    if( tools.contains("strvctvre") || tools.contains("strvctvre") ) {
        STRVCTVRE_STRVCTVRE(
            ch_svannotate.map { meta, vcf, tbi -> [ meta, vcf, tbi, genome ] },
            strvctvre_phylop,
            strvctvre_data
        )
        ch_versions = ch_versions.mix(STRVCTVRE_STRVCTVRE.out.versions.first())

        TABIX_BGZIPTABIX(
            STRVCTVRE_STRVCTVRE.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

        ch_strvctvre = TABIX_BGZIPTABIX.out.gz_tbi
    }

    emit:
    vcfs = ch_strvctvre
    versions = ch_versions
    reports = ch_reports
}
