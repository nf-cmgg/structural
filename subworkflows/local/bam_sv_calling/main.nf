//
// Call structural variants
//

// Import subworkflows
include { BAM_VARIANT_CALLING_MANTA                     } from '../bam_variant_calling_manta/main'
include { BAM_VARIANT_CALLING_DELLY                     } from '../bam_variant_calling_delly/main'
include { BAM_VARIANT_CALLING_WHAMG                     } from '../bam_variant_calling_whamg/main'
include { BAM_VARIANT_CALLING_SMOOVE                    } from '../bam_variant_calling_smoove/main'
include { BAM_VARIANT_CALLING_SCRAMBLE                  } from '../bam_variant_calling_scramble/main'
include { BAM_VARIANT_CALLING_GRIDSS                    } from '../bam_variant_calling_gridss/main'
include { VCF_MERGE_JASMINE                             } from '../vcf_merge_jasmine/main'

// Import modules
include { REHEADER_CALLED_VCFS                          } from '../../../modules/local/bcftools/reheader_called_vcfs/main'

include { BCFTOOLS_SORT                                 } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX                                   } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_SV_CALLING {
    take:
        ch_crams        // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta        // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai          // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_bwa_index    // channel: [optional]  [ meta, index ] => The BWA MEM index

    main:

    val_callers     = params.callers.tokenize(",").intersect(params.sv_callers)

    ch_versions     = Channel.empty()
    ch_reports      = Channel.empty()
    ch_called_vcfs  = Channel.empty()

    //
    // Calling variants using Manta
    //

    if("manta" in val_callers){
        BAM_VARIANT_CALLING_MANTA(
            ch_crams,
            ch_fasta,
            ch_fai
        )

        ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_MANTA.out.manta_vcfs)
        ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_MANTA.out.versions)
    }

    //
    // Calling variants using Delly
    //

    if("delly" in val_callers){
        BAM_VARIANT_CALLING_DELLY(
            ch_crams,
            ch_fasta,
            ch_fai
        )

        ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_DELLY.out.delly_vcfs)
        ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_DELLY.out.versions)
    }

    //
    // Calling variants using Whamg (Currently disabled)
    //

    // TODO Whamg needs some reheadering (like done in https://github.com/broadinstitute/gatk-sv/blob/90e3e9a221bdfe7ab2cfedeffb704bc6f0e99aa9/wdl/Whamg.wdl#L209)
    // TODO Add insertions sequence in the info key - Whamg will not work for now
    // if("whamg" in val_callers){
    //     BAM_VARIANT_CALLING_WHAMG(
    //         ch_crams,
    //         ch_fasta,
    //         ch_fai
    //     )

    //     ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_WHAMG.out.whamg_vcfs)
    //     ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_WHAMG.out.versions)
    // }

    //
    // Calling variants using Smoove
    //

    if("smoove" in val_callers){
        BAM_VARIANT_CALLING_SMOOVE(
            ch_crams,
            ch_fasta,
            ch_fai
        )

        ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_SMOOVE.out.smoove_vcfs)
        ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_SMOOVE.out.versions)
    }

    //
    // Calling variants using Gridss
    //

    if("gridss" in val_callers){
        BAM_VARIANT_CALLING_GRIDSS(
            ch_crams,
            ch_fasta,
            ch_fai,
            ch_bwa_index
        )

        ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_GRIDSS.out.gridss_vcfs)
        ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_GRIDSS.out.versions)
    }

    //
    // Calling variants using Scramble (I don't know if calling variants is the correct term here)
    //

    // Scramble is unfinished. It needs a lot of improvements if we were to add it

    // if("scramble" in val_callers){
    //     BAM_VARIANT_CALLING_SCRAMBLE(
    //         ch_crams,
    //         ch_fasta
    //     )

    //    ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_SCRAMBLE.out.scramble_vcfs)
    //    ch_versions     = ch_versions.mix(BAM_VARIANT_CALLING_SCRAMBLE.out.versions)
    // }

    //
    // Standardize and merge VCFs per sample for all callers
    //

    ch_called_vcfs
        .map { meta, vcf, tbi ->
            [ meta, vcf ]
        }
        .set { ch_merge_input }

    VCF_MERGE_JASMINE(
        ch_merge_input,
        ch_fasta,
        ch_fai,
    )
    ch_versions = ch_versions.mix(VCF_MERGE_JASMINE.out.versions)

    VCF_MERGE_JASMINE.out.merged_vcfs.set { ch_merged_vcfs }

    emit:
    vcfs                = ch_merged_vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions            = ch_versions
    reports             = ch_reports
}