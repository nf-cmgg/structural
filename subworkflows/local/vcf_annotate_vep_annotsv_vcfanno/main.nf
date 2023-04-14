//
// Annotate the VCFs
//

include { ANNOTSV_ANNOTSV                       } from '../../../modules/nf-core/annotsv/annotsv/main'
include { ENSEMBLVEP_VEP                        } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                               } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIPTABIX as TABIX_ANNOTATED   } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VEP              } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO {
    take:
        ch_vcfs                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
        ch_fasta                // channel: [mandatory] [ path(fasta) ] => The fasta reference file
        ch_fai                  // channel: [mandatory] [ path(fai) ] => The index of the fasta reference file
        ch_annotsv_annotations  // channel: [mandatory] [ val(meta), path(annotations) ] => The annotations for AnnotSV
        ch_vep_cache            // channel: [optional]  [ path(cache) ] => The path to the local VEP cache
        ch_vep_extra_files      // channel: [optional]  [ path(file1, file2, file3...) ] => The VEP extra files
        ch_vcfanno_toml         // channel: [mandatory] [ path(toml) ] => The TOML configuration for VCFanno
        ch_vcfanno_lua          // channel: [optional]  [ path(lua) ] => The lua script to influence VCFanno
        ch_vcfanno_resources    // channel: [optional]  [ path(file1, file2, file3...) ] => The extra VCFanno files

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // Run AnnotSV and VEP in parallel and merge TSV from AnnotSV with VCF from VEP during VCFanno

    // TODO add the other inputs
    ANNOTSV_ANNOTSV(
        ch_vcfs,
        ch_annotsv_annotations,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions)

    ENSEMBLVEP_VEP(
        ch_vcfs,
        params.genome,
        params.species,
        params.vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        ch_vep_extra_files
    )
    ch_reports  = ch_reports.mix(ENSEMBLVEP_VEP.out.report)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)

    ANNOTSV_ANNOTSV.out.tsv
        .join(ENSEMBLVEP_VEP.out.vcf, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_vcfanno_input }

    // TODO update vcfanno to work with sample specific custom files
    VCFANNO(
        ch_vcfanno_input,
        ch_vcfanno_toml,
        ch_vcfanno_lua,
        ch_vcfanno_resources
    )
    ch_versions = ch_versions.mix(VCFANNO.out.versions)

    TABIX_ANNOTATED(
        VCFANNO.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_ANNOTATED.out.versions)
    
    emit:
    annotated_vcfs  = TABIX_ANNOTATED.out.gz_tbi  // channel: [ val(meta), path(vcf), path(tbi) ]

    reports         = ch_reports
    versions        = ch_versions
}
