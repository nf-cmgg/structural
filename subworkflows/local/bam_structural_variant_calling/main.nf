//
// Gather Sample Evidence
//

// Import subworkflows
include { BAM_VARIANT_CALLING_MANTA                     } from '../bam_variant_calling_manta/main'
include { BAM_VARIANT_CALLING_DELLY                     } from '../bam_variant_calling_delly/main'
include { BAM_VARIANT_CALLING_WHAMG                     } from '../bam_variant_calling_whamg/main'
include { BAM_VARIANT_CALLING_SMOOVE                    } from '../bam_variant_calling_smoove/main'
include { BAM_VARIANT_CALLING_SCRAMBLE                  } from '../bam_variant_calling_scramble/main'
include { VCF_METRICS_SVTK_SVTEST                       } from '../vcf_metrics_svtk_svtest/main'
include { VCF_MERGE_JASMINE                             } from '../vcf_merge_jasmine/main'

include { BAM_VARIANT_CALLING_GRIDSS                    } from '../bam_variant_calling_gridss/main'

// Import modules
include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS  } from '../../../modules/nf-core/gatk4/collectreadcounts/main'
include { GATK4_COLLECTSVEVIDENCE as COLLECTSVEVIDENCE  } from '../../../modules/nf-core/gatk4/collectsvevidence/main'
include { TABIX_TABIX                                   } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_STRUCTURAL_VARIANT_CALLING {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        allele_loci_vcf         // channel: [optional]  [ vcf ] => A channel containing the VCF and its index for counting the alleles
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file
        dict                    // channel: [mandatory] [ dict ] => The dictionary of the fasta reference file
        bwa_index               // channel: [optional]  [ index ] => The BWA MEM index

    main:

    callers = params.callers.tokenize(",")

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    ch_metrics  = Channel.empty()
    called_vcfs = Channel.empty()

    //
    // GATK Collect Read Counts
    //

    crams
        .combine(
            beds.map({meta, bed, bed_gz, bed_gz_tbi -> [meta, bed]})
        , by:0)
    .dump(tag: 'collectreadcounts_input', pretty: true)
    .set { collectreadcounts_input }

    COLLECTREADCOUNTS(
        collectreadcounts_input,
        fasta,
        fasta_fai,
        dict
    )

    // ch_versions = ch_versions.mix(COLLECTREADCOUNTS.out.versions)

    //
    // Calling variants using Manta
    //

    if("manta" in callers){
        BAM_VARIANT_CALLING_MANTA(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(BAM_VARIANT_CALLING_MANTA.out.manta_vcfs)
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_MANTA.out.versions)
    }

    //
    // Calling variants using Delly
    //

    if("delly" in callers){
        BAM_VARIANT_CALLING_DELLY(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(BAM_VARIANT_CALLING_DELLY.out.delly_vcfs)
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_DELLY.out.versions)
    }

    //
    // Calling variants using Whamg
    //

    // Whamg needs some reheadering (like done in https://github.com/broadinstitute/gatk-sv/blob/90e3e9a221bdfe7ab2cfedeffb704bc6f0e99aa9/wdl/Whamg.wdl#L209)
    if("whamg" in callers){
        BAM_VARIANT_CALLING_WHAMG(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(BAM_VARIANT_CALLING_WHAMG.out.whamg_vcfs)
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_WHAMG.out.versions)
    }

    //
    // Calling variants using Smoove
    //

    if("smoove" in callers){
        BAM_VARIANT_CALLING_SMOOVE(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(BAM_VARIANT_CALLING_SMOOVE.out.smoove_vcfs)
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SMOOVE.out.versions)
    }

    //
    // Calling variants using Gridss
    //

    if("gridss" in callers){
        BAM_VARIANT_CALLING_GRIDSS(
            crams,
            fasta,
            fasta_fai,
            bwa_index
        )

        called_vcfs = called_vcfs.mix(BAM_VARIANT_CALLING_GRIDSS.out.gridss_vcfs)
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_GRIDSS.out.versions)
    }

    //
    // Calling variants using Scramble (I don't know if calling variants is the correct term here)
    //

    // Scramble is unfinished. It needs a lot of improvements if we were to add it

    // if("scramble" in callers){
    //     BAM_VARIANT_CALLING_SCRAMBLE(
    //         crams,
    //         beds,
    //         fasta
    //     )

    //     called_vcfs = called_vcfs.mix(BAM_VARIANT_CALLING_SCRAMBLE.out.scramble_vcfs)
    //     ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_SCRAMBLE.out.versions)
    // }

    //
    // GATK Collect Structural Variant Evidence
    //

    if(allele_loci_vcf){
        TABIX_TABIX(
            allele_loci_vcf.map{[[id:"allele_loci_vcf"], it]}
        )

        crams
            .combine(TABIX_TABIX.out.tbi)
            .map(
                { meta, cram, crai, tbi ->
                    [ meta, cram, crai, allele_loci_vcf, tbi ]
                }
            )
            .set { collectsvevidence_input }

    } else {
        crams
            .map(
                { meta, cram, crai ->
                    [ meta, cram, crai, [], [] ]
                }
            )
            .set { collectsvevidence_input }
    }

    collectsvevidence_input.dump(tag: 'collectsvevidence_input', pretty: true)

    COLLECTSVEVIDENCE(
        collectsvevidence_input,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(COLLECTSVEVIDENCE.out.versions)

    //
    // Create the metrics for all produced files
    //

    if(params.run_module_metrics) {
        VCF_METRICS_SVTK_SVTEST(
            called_vcfs,
            COLLECTSVEVIDENCE.out.split_read_evidence,
            COLLECTSVEVIDENCE.out.paired_end_evidence,
            COLLECTSVEVIDENCE.out.site_depths,
            fasta_fai
        )

        ch_metrics  = ch_reports.mix(VCF_METRICS_SVTK_SVTEST.out.metrics)
        ch_versions = ch_versions.mix(VCF_METRICS_SVTK_SVTEST.out.versions)
    }

    if(callers.size > 1){
        VCF_MERGE_JASMINE(
            called_vcfs.map { it[0..1] },
            fasta,
            fasta_fai,
        )

        VCF_MERGE_JASMINE.out.merged_vcfs.set { merged_vcfs }
    } else {
        called_vcfs.set { merged_vcfs }
    }

    emit:
    vcfs                = called_vcfs

    coverage_counts     = COLLECTREADCOUNTS.out.tsv

    split_reads         = COLLECTSVEVIDENCE.out.split_read_evidence.combine(COLLECTSVEVIDENCE.out.split_read_evidence_index, by:0)
    read_pairs          = COLLECTSVEVIDENCE.out.paired_end_evidence.combine(COLLECTSVEVIDENCE.out.paired_end_evidence_index, by:0)
    site_depths         = COLLECTSVEVIDENCE.out.site_depths.combine(COLLECTSVEVIDENCE.out.site_depths_index, by:0)

    versions            = ch_versions
    metrics             = ch_metrics
    reports             = ch_reports
}
