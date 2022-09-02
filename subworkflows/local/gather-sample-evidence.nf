//
// Gather Sample Evidence
//
include { WHAMG            } from '../../modules/nf-core/modules/whamg/main'
include { DELLY_CALL       } from '../../modules/nf-core/modules/delly/call/main'
include { BCFTOOLS_CONVERT } from '../../modules/nf-core/modules/bcftools/convert/main'
include { MANTA_GERMLINE   } from '../../modules/nf-core/modules/manta/germline/main'
include { SAMTOOLS_CONVERT } from '../../modules/nf-core/modules/samtools/convert/main'

workflow GATHER_SAMPLE_EVIDENCE {
    take:
        crams       // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        beds        // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta       // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai   // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file
        dict        // channel: [mandatory] [ dict ] => The dictionary of the fasta reference file

        callers     // array:   [mandatory] => A list of callers to use for variant calling, should be one of these: delly|manta|whamg

    main:

    called_vcfs = Channel.empty()

    if("delly" in callers){
        DELLY_CALL(
            crams,
            fasta,
            fasta_fai
        )

        BCFTOOLS_CONVERT(
            DELLY_CALL.out.bcf.map({ meta, bcf -> [ meta, bcf, [] ]}),
            [],
            fasta
        )

        called_vcfs = called_vcfs.mix(
            BCFTOOLS_CONVERT.out.vcf_gz.map({ meta, vcf -> 
                new_meta = meta.clone()
                new_meta.caller = "delly"
                [ new_meta, vcf ]
            })
        )

    }

    if("manta" in callers){
        // TODO add a check for germline or somatic
        manta_input = crams.combine(beds, by: 0).map({ meta, cram, crai, bed, bed_gz, bed_gz_tbi ->
            [ meta, cram, crai, bed_gz, bed_gz_tbi ]
        })

        MANTA_GERMLINE(
            manta_input,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(
            MANTA_GERMLINE.out.diploid_sv_vcf.map({ meta, vcf -> 
                new_meta = meta.clone()
                new_meta.caller = "manta"
                [ new_meta, vcf ]
            })
        )        
    }

    if("whamg" in callers){
        SAMTOOLS_CONVERT(
            crams,
            fasta,
            fasta_fai
        )

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
    }

    emit:
    vcfs = called_vcfs

}