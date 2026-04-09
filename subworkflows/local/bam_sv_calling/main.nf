//
// Call structural variants
//

// Import subworkflows
include { BAM_VARIANT_CALLING_MANTA  } from '../bam_variant_calling_manta/main'
include { BAM_VARIANT_CALLING_DELLY  } from '../bam_variant_calling_delly/main'
include { BAM_VARIANT_CALLING_SMOOVE } from '../bam_variant_calling_smoove/main'
include { VCF_MERGE_CALLERS_JASMINE  } from '../vcf_merge_callers_jasmine/main'


workflow BAM_SV_CALLING {
    take:
        ch_crams            // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta            // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai              // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_manta_config     // channel: [optional]  [ config ] => The config to pass to Manta
        ch_svync_configs    // channel: [mandatory] [ configs ] => A list of svync config files
        val_callers         // value:   [mandatory] => List of all SV callers to use

    main:
    def ch_reports      = channel.empty()
    def ch_called_vcfs  = channel.empty()
    def ch_raw_vcfs     = channel.empty()

    //
    // Calling variants using Manta
    //

    if("manta" in val_callers){
        BAM_VARIANT_CALLING_MANTA(
            ch_crams,
            ch_fasta,
            ch_fai,
            ch_manta_config,
            ch_svync_configs
        )

        ch_raw_vcfs     = ch_raw_vcfs.mix(BAM_VARIANT_CALLING_MANTA.out.raw_vcfs)
        ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_MANTA.out.manta_vcfs)
    }

    //
    // Calling variants using Delly
    //

    if("delly" in val_callers){
        BAM_VARIANT_CALLING_DELLY(
            ch_crams,
            ch_fasta,
            ch_fai,
            ch_svync_configs
        )

        ch_raw_vcfs     = ch_raw_vcfs.mix(BAM_VARIANT_CALLING_DELLY.out.raw_vcfs)
        ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_DELLY.out.delly_vcfs)
    }

    //
    // Calling variants using Smoove
    //

    if("smoove" in val_callers){
        BAM_VARIANT_CALLING_SMOOVE(
            ch_crams,
            ch_fasta,
            ch_fai,
            ch_svync_configs
        )

        ch_raw_vcfs     = ch_raw_vcfs.mix(BAM_VARIANT_CALLING_SMOOVE.out.raw_vcfs)
        ch_called_vcfs  = ch_called_vcfs.mix(BAM_VARIANT_CALLING_SMOOVE.out.smoove_vcfs)
    }

    def ch_merged_vcfs = channel.empty()
    if(val_callers.size() > 1) {
        VCF_MERGE_CALLERS_JASMINE(
            ch_called_vcfs,
            ch_fasta,
            ch_fai,
            val_callers,
            "sv"
        )
        ch_merged_vcfs = VCF_MERGE_CALLERS_JASMINE.out.vcfs
    } else {
        ch_merged_vcfs = ch_called_vcfs
            .map { meta, vcf, tbi ->
                def new_meta = meta - meta.subMap("caller") + [variant_type:"sv"]
                [ new_meta, vcf, tbi ]
            }
    }

    emit:
    caller_vcfs         = ch_raw_vcfs       // channel: [ val(meta), path(vcf), path(tbi)]
    vcfs                = ch_merged_vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]

    reports             = ch_reports
}
