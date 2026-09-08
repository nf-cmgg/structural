//
// Call variants in BAM files using SMN Copy Number Caller
//

include { SMNCOPYNUMBERCALLER   } from '../../../modules/nf-core/smncopynumbercaller/main'

workflow BAM_VARIANT_CALLING_SMNCOPYNUMBERCALLER {
    take:
        ch_crams        // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta        // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai          // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file

    main:
    def ch_smncopynumbercaller_input = ch_crams
        .map { meta, cram, crai ->
            def new_meta = meta + [ caller:'smncopynumbercaller' ]
            [ new_meta, cram, crai ]
        }

    SMNCOPYNUMBERCALLER(
        ch_smncopynumbercaller_input,
        ch_fasta.join(ch_fai, failOnDuplicate:true, failOnMismatch:true).collect()
    )

    def ch_smncopynumbercaller_out = SMNCOPYNUMBERCALLER.out.smncopynumber
            .mix(SMNCOPYNUMBERCALLER.out.run_metrics)

    emit:
    caller_out = ch_smncopynumbercaller_out   // channel: [ val(meta), path(file) ] => The output of the SMN Copy Number Caller
}
