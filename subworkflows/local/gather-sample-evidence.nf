//
// Gather Sample Evidence
//

// Import subworkflows
include { RUN_DELLY                                     } from './run-delly'
include { RUN_MANTA                                     } from './run-manta'

// Import modules
include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS  } from '../../modules/nf-core/modules/gatk4/collectreadcounts/main'
include { WHAMG                                         } from '../../modules/nf-core/modules/whamg/main'
include { MANTA_GERMLINE                                } from '../../modules/nf-core/modules/manta/germline/main'
include { SAMTOOLS_CONVERT                              } from '../../modules/nf-core/modules/samtools/convert/main'

workflow GATHER_SAMPLE_EVIDENCE {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file
        dict                    // channel: [mandatory] [ dict ] => The dictionary of the fasta reference file

    main:

    callers = params.callers.split(",")

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    called_vcfs = Channel.empty()

    collectreadcounts_input = crams.combine(
        beds.map({meta, bed, bed_gz, bed_gz_tbi -> [meta, bed]}), by:0
    )

    COLLECTREADCOUNTS(
        collectreadcounts_input,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(COLLECTREADCOUNTS.out.versions)

    if("delly" in callers){
        RUN_DELLY(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(RUN_DELLY.out.delly_vcfs)
        ch_versions = ch_versions.mix(RUN_DELLY.out.versions)
    }

    if("manta" in callers){
        RUN_MANTA(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(RUN_MANTA.out.manta_vcfs)
        ch_versions = ch_versions.mix(RUN_MANTA.out.versions)
    }

    if("whamg" in callers){
        SAMTOOLS_CONVERT(
            crams,
            fasta,
            fasta_fai
        )
        ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

        WHAMG(
            SAMTOOLS_CONVERT.out.alignment_index,
            fasta
        )

        called_vcfs = called_vcfs.mix(
            WHAMG.out.vcf.map({ meta, vcf -> 
                new_meta = meta.clone()
                new_meta.caller = "whamg"
                [ new_meta, vcf ]
            })
        )        
        ch_versions = ch_versions.mix(WHAMG.out.versions)
    }

    emit:
    vcfs            = called_vcfs.view()
    coverage_counts = COLLECTREADCOUNTS.out.tsv
    versions        = ch_versions
    reports         = ch_reports
}