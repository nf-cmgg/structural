//
// Gather Batch Evidence
//

// Import subworkflows
include { MAKE_BINCOV_MATRIX                            } from '../common-workflows/make-bincov-matrix'
include { BATCH_EVIDENCE_MERGING                        } from './batch-evidence-merging'

// Import modules
include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS  } from '../../../modules/nf-core/gatk4/collectreadcounts/main'

workflow GATHER_BATCH_EVIDENCE {
    take:
        count_files             // channel: [mandatory] [ meta, count_file ] => The read counts of all samples
        bincov_matrix           // channel: [optional]  [ bincov_matrix ] => The bin coverage matrix
        bincov_matrix_index     // channel: [optional]  [ bincov_matrix_index ] => The bin coverage matrix index

        BAF_files               // channel: [optional]  [ baf_files, index ] => The BAF files
        PE_files                // channel: [mandatory] [ meta, pe_file, index ] => The paired end evidence files created in GatherSampleEvidence
        SR_files                // channel: [mandatory] [ meta, sr_file, index ] => The split read evidence files created in GatherSampleEvidence

        SD_files                // channel: [optional]  [ meta, sd_file, index ] => The site depth evidence files created in GatherSampleEvidence
        allele_loci_vcf         // channel: [optional]  [ vcf, tbi ] => VCF of SNPs marking loci for allele count

        fasta                   // channel: [mandatory] [ fasta ] => The reference FASTA file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the FASTA reference file
        dict                    // channel: [mandatory] [ dict ] => The sequence dictionary of the reference genome


    main:

    ch_versions         = Channel.empty()
    batch_ploidy_matrix = Channel.empty()
    batch_ploidy_plots  = Channel.empty()

    //
    // Bin coverage matrix
    //

    if(!bincov_matrix || !bincov_matrix_index){

        MAKE_BINCOV_MATRIX(
            count_files
        )

        ch_versions = ch_versions.mix(MAKE_BINCOV_MATRIX.out.versions)

        MAKE_BINCOV_MATRIX.out.merged_bincov.set { bincov_matrix }
        MAKE_BINCOV_MATRIX.out.merged_bincov_index.set { bincov_matrix_index }

    }

    //
    // Ploidy
    //

    if(params.run_ploidy){

        PLOIDY(
            bincov_matrix
        )

        ploidy_matrix   = ploidy_matrix.mix(PLOIDY.out.ploidy_matrix)
        ploidy_plots    = ploidy_plots.mix(PLOIDY.out.ploidy_plots)
        ch_versions     = ch_versions.mix(PLOIDY.out.versions)

    }

    // SubsetPedFile (this will be left out for now -> subsets the PED file so that it only contains all present samples)

    // AddCaseSampleToBed (optional => adds the first sample to PED)

    //
    // BatchEvidenceMerging
    //

    BATCH_EVIDENCE_MERGING(
        BAF_files,
        PE_files,
        SR_files,
        SD_files,
        allele_loci_vcf,
        fasta,
        fasta_fai,
        dict
    )

    // CNMOPS (CNV related)

    // CNMOPSLarge (CNV related)

    // CondenseReadCounts (uses GATK CondenseDepthEvidence, which is not in a released version yet)

    // gCNVCase (CNV related)

    // MergeDepth (CNV related)

    // MedianCov

    // PreprocessPESR

    // TinyResolve

    // MatrixQC

    // GatherBatchEvidenceMetrics

    emit:
    versions            = ch_versions
}