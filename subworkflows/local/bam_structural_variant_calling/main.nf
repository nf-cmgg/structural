//
// Gather Sample Evidence
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
include { VIOLA                                         } from '../../../modules/local/viola/main'
include { REHEADER_CALLED_VCFS                          } from '../../../modules/local/bcftools/reheader_called_vcfs/main'

include { TABIX_TABIX as TABIX_VCFS                     } from '../../../modules/nf-core/tabix/tabix/main'

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
    // Calling variants using Whamg (Currently disabled)
    //

    // TODO Whamg needs some reheadering (like done in https://github.com/broadinstitute/gatk-sv/blob/90e3e9a221bdfe7ab2cfedeffb704bc6f0e99aa9/wdl/Whamg.wdl#L209)
    // TODO Add insertions sequence in the info key - Whamg will not work for now
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
    // Calling variants using Gridss (Currently disabled)
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
    // Standardize and merge VCFs per sample for all callers
    //

    VIOLA(
        called_vcfs.map{ it[0..1] }
    )

    ch_versions = ch_versions.mix(VIOLA.out.versions)

    if(callers.size > 1){
        VCF_MERGE_JASMINE(
            VIOLA.out.vcf,
            fasta,
            fasta_fai,
        )
        ch_versions = ch_versions.mix(VCF_MERGE_JASMINE.out.versions)

        VCF_MERGE_JASMINE.out.merged_vcfs.set { merged_vcfs }
    } else {

        new_header = Channel.fromPath("${projectDir}/assets/header.txt").collect()

        REHEADER_CALLED_VCFS(
            VIOLA.out.vcf,
            new_header,
            fasta_fai
        )
        ch_versions = ch_versions.mix(REHEADER_CALLED_VCFS.out.versions)
        
        TABIX_VCFS(
            REHEADER_CALLED_VCFS.out.vcf
        )
        ch_versions = ch_versions.mix(TABIX_VCFS.out.versions)

        REHEADER_CALLED_VCFS.out.vcf
            .join(TABIX_VCFS.out.tbi, failOnDuplicate:true, failOnMismatch:true)
            .map { meta, vcf, tbi ->
                new_meta = meta - meta.subMap("caller")
                [ new_meta, vcf, tbi ]
            }
            .set { merged_vcfs }
    }

    emit:
    vcfs                = merged_vcfs

    versions            = ch_versions
    metrics             = ch_metrics
    reports             = ch_reports
}
