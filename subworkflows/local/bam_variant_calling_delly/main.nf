nextflow.preview.types = true
//
// Run Delly
//

include { DELLY_CALL; DellyCallInput  } from '../../../modules/nf-core/delly/call/main'
include { BCFTOOLS_CONCAT; BcftoolsConcatInput } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT     } from '../../../modules/nf-core/bcftools/sort/main'
include { SVYNC; SvyncInput } from '../../../modules/nf-core/svync/main'

workflow BAM_VARIANT_CALLING_DELLY {
    take:
    ch_input: Channel<DellyCallInput> // The aligned CRAMs per sample with the regions they should be called on

    main:

    //
    // Calling variants using Delly
    //

    def ch_delly_input: Channel<DellyCallInput> = ch_input
        .combine(channel.of("DEL", "INS", "INV", "DUP", "BND"))
        .map { rec, sv_type ->
            rec + record(sv_type: sv_type)
        }

    def delly_out = DELLY_CALL(
        ch_delly_input
    )

    def ch_concat_input: Channel<BcftoolsConcatInput> = delly_out
        .map { rec: Record ->
            tuple(rec.id, 5, rec)
        }
        .groupBy()
        .map { _id, recs ->
            // Workaround for TaskPath issue: https://github.com/nextflow-io/nextflow/issues/7032
            def vcfs = []
            def tbis = []
            recs.each { r ->
                vcfs << r.bcf.toRealPath()
                tbis << r.csi.toRealPath()
            }
            def first_rec = recs.toList().first()
            first_rec + record(
                vcfs: vcfs,
                tbis: tbis,
                input: first_rec.input.toRealPath(),
                input_index: first_rec.input_index.toRealPath(),
                fasta: first_rec.fasta.toRealPath(),
                fai: first_rec.fai.toRealPath(),
            )
        }

    def concat_out = BCFTOOLS_CONCAT(
        ch_concat_input
    )

    def sort_out = BCFTOOLS_SORT(
        concat_out.map { rec ->
            rec + record(
                vcfs: rec.vcfs,
                tbis: rec.tbis,
            )
        }
    )

    def ch_svync_input: Channel<SvyncInput> = sort_out
        .combine(config:file("${projectDir}/assets/svync/delly.yaml"))

    def svync_out = SVYNC(
        ch_svync_input
    )

    emit:
    raw_vcfs    = sort_out // channel: [ val(meta), path(vcf), path(tbi) ]
    delly_vcfs  = svync_out   // channel: [ val(meta), path(vcf), path(tbi) ]
}
