/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_structural_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BAM_PREPARE_SAMTOOLS                  } from '../subworkflows/local/bam_prepare_samtools/main'
include { BAM_SV_CALLING                        } from '../subworkflows/local/bam_sv_calling/main'
include { BAM_CNV_CALLING                       } from '../subworkflows/local/bam_cnv_calling/main'
include { VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO      } from '../subworkflows/local/vcf_annotate_vep_annotsv_vcfanno/main'
include { BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER } from '../subworkflows/local/bam_repeat_estimation_expansionhunter/main'
include { VCF_CONCAT_BCFTOOLS                   } from '../subworkflows/local/vcf_concat_bcftools/main'
include { VCF_MERGE_FAMILY_JASMINE              } from '../subworkflows/local/vcf_merge_family_jasmine/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_FAIDX                    } from '../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                         } from '../modules/nf-core/bwa/index/main'
include { ENSEMBLVEP_DOWNLOAD               } from '../modules/nf-core/ensemblvep/download/main'
include { ANNOTSV_INSTALLANNOTATIONS        } from '../modules/nf-core/annotsv/installannotations/main'
include { UNTAR as UNTAR_ANNOTSV            } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BWA                } from '../modules/nf-core/untar/main'
include { NGSBITS_SAMPLEGENDER              } from '../modules/nf-core/ngsbits/samplegender/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STRUCTURAL {

    take:
    // channels
    ch_samplesheet              // channel: samplesheet read in from --input
    ch_multiqc_config
    ch_multiqc_custom_config
    ch_multiqc_logo
    ch_multiqc_custom_methods_description

    // file inputs
    fasta                       // The fasta file to use
    fai                         // The index of the fasta file
    expansionhunter_catalog     // The expansionhunter catalog
    qdnaseq_female              // The QDNAseq annotations for female samples
    qdnaseq_male                // The QDNAseq annotations for male samples
    wisecondorx_reference       // The WisecondorX annotations file
    vep_cache                   // The VEP cache directory
    annotsv_annotations         // The annotations directory for AnnotSV
    annotsv_candidate_genes     // A file containing the AnnotSV candidate genes
    annotsv_gene_transcripts    // A file containing the AnnotSV gene transcripts
    vcfanno_lua                 // A Lua script to use with vcfanno
    vcfanno_resources           // A comma delimited list of paths to vcfanno resource files
    vcfanno_toml                // The vcfanno config to
    blacklist                   // A BED file of blacklisted regions
    manta_config                // A configuration file to be used in Manta
    svync_dir                   // A directory containing svync configs (they must end in '.yaml' and contain the full name of the caller)
    bedgovcf_dir                // A directory containing bedgovcf configs (they must end in '.yaml' and contain the full name of the caller)
    vcfanno_default_dir         // A directory containing the default vcfanno configs (they must end in '.toml')

    // boolean inputs
    annotate                    // Run annotation on SV and CNV data
    concat_output               // Concatenate all output files

    // value inputs
    input_callers               // All callers to be used
    sv_callers_support          // How many SV callers are needed to support the variant
    cnv_callers_support         // How many CNV callers are needed to support the variant
    genome                      // The genome to use
    species                     // The species to be used by VEP
    vep_assembly                // The genome assembly to be downloaded for VEP
    vep_cache_version           // The version of the VEP cache to use
    annotations_filter          // The filter pattern to use after annotation
    outdir                      // The output directory of the pipeline

    main:

    ch_versions         = Channel.empty()
    ch_reports          = Channel.empty()
    ch_outputs          = Channel.empty()
    ch_annotation_input = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    variant_types = [] // The variant types that can be annotated this run
    count_types = 0 // The amount of different variant types that can be concatenated

    //
    // Input validation
    //

    // When making changes here, make sure to also update the following files: conf/modules.config
    def svCallers = ["delly", "manta", "smoove"] //, "gridss"
    def cnvCallers = ["qdnaseq", "wisecondorx"]
    def repeatsCallers = ["expansionhunter"]

    def allCallers = svCallers + cnvCallers + repeatsCallers
    def annotationCallers = svCallers + cnvCallers

    // Callers that need the sex
    def sexCallers = ["expansionhunter", "qdnaseq"]

    def callers = input_callers.toLowerCase().tokenize(",").collect { caller ->
        if(caller == "all") {return allCallers}
        if(caller == "sv")  {return svCallers}
        if(caller == "cnv") {return cnvCallers}
        if(caller == "rre") {return repeatsCallers}
        return caller
    }.flatten()

    callers.each { caller ->
        if(!(caller in allCallers)) { error("The caller '${caller}' is not supported please specify a comma delimited list with on or more of the following callers: ${allCallers}".toString()) }
    }

    def sv_callers_to_use = callers.intersect(svCallers)
    def cnv_callers_to_use = callers.intersect(cnvCallers)
    def rre_callers_to_use = callers.intersect(repeatsCallers)

    if (sv_callers_to_use && sv_callers_support > sv_callers_to_use.size()) {
        error("The 'sv_callers_support' option (${sv_callers_support}) is higher than the amount of SV callers given (${sv_callers_to_use.size()}). Please adjust 'sv_callers_support' to a value lower of equal to the amount of SV callers to use.")
    }

    if (cnv_callers_to_use && cnv_callers_support > cnv_callers_to_use.size()) {
        error("The 'cnv_callers_support' option (${cnv_callers_support}) is higher than the amount of CNV callers given (${cnv_callers_to_use.size()}). Please adjust 'cnv_callers_support' to a value lower of equal to the amount of CNV callers to use.")
    }

    if ("qdnaseq" in callers && (!qdnaseq_male || !qdnaseq_female)) {
        error("Please give the QDNAseq references using 'qdnaseq_male' and 'qdnaseq_female'")
    }

    if ("wisecondorx" in callers && !wisecondorx_reference) {
        error("Please give the WisecondorX reference using 'wisecondorx_reference'")
    }

    //
    // Create input channels from parameters
    //

    ch_fasta                    = Channel.fromPath(fasta).collect { fasta_file -> [[id:'fasta'], fasta_file ] }
    ch_annotsv_candidate_genes  = annotsv_candidate_genes ?  Channel.fromPath(annotsv_candidate_genes).collect { genes_file -> [[], genes_file] } : [[],[]]
    ch_annotsv_gene_transcripts = annotsv_gene_transcripts ? Channel.fromPath(annotsv_gene_transcripts).collect { transcripts_file -> [[], transcripts_file] } : [[],[]]
    ch_vcfanno_lua              = vcfanno_lua ?              Channel.fromPath(vcfanno_lua).collect() : []
    ch_catalog                  = expansionhunter_catalog ?  Channel.fromPath(expansionhunter_catalog).collect { catalog_file -> [[id:'catalog'], catalog_file] } : [[],[]]
    ch_qdnaseq_male             = qdnaseq_male ?             Channel.fromPath(qdnaseq_male).collect { qdnaseq_file -> [[id:'qdnaseq_male'], qdnaseq_file] } : [[],[]]
    ch_qdnaseq_female           = qdnaseq_female ?           Channel.fromPath(qdnaseq_female).collect { qdnaseq_file -> [[id:'qdnaseq_female'], qdnaseq_file] } : [[],[]]
    ch_wisecondorx_reference    = wisecondorx_reference ?    Channel.fromPath(wisecondorx_reference).collect { wcx_file -> [[id:'wisecondorx'], wcx_file] } : [[],[]]
    ch_blacklist                = blacklist ?                Channel.fromPath(blacklist).collect { blacklist_file -> [[id:'blacklist'], blacklist_file] } : [[],[]]
    ch_manta_config             = manta_config ?             Channel.fromPath(manta_config).collect() : [[]]
    ch_svync_configs            = svync_dir ?                Channel.fromPath("${svync_dir}/*.yaml", checkIfExists:true).collect() : []
    ch_bedgovcf_configs         = bedgovcf_dir ?             Channel.fromPath("${bedgovcf_dir}/*.yaml", checkIfExists:true).collect() : []

    val_vcfanno_resources       = vcfanno_resources ?        vcfanno_resources.split(",").collect { resource_file -> file(resource_file, checkIfExists:true) }.flatten() : []
    val_default_vcfanno_tomls   = vcfanno_default_dir ?      files("${vcfanno_default_dir}/*.toml", checkIfExists:true) : []
    val_vcfanno_toml            = vcfanno_toml ?             file(vcfanno_toml, checkIfExists:true) : []

    ch_vep_extra_files = []

    //
    // Create optional inputs
    //

    if(!fai){
        SAMTOOLS_FAIDX(
            ch_fasta
        )

        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai      = SAMTOOLS_FAIDX.out.fai.collect { fai_file -> [[id:'fai'], fai_file] }
    }
    else {
        ch_fai = Channel.fromPath(fai).collect { fai_file -> [[id:'fai'], fai_file] }
    }

    // if(!bwa && "gridss" in callers){
    //     BWA_INDEX(
    //         ch_fasta
    //     )

    //     ch_versions  = ch_versions.mix(BWA_INDEX.out.versions)
    //     ch_bwa_index = BWA_INDEX.out.index.map{[[id:'bwa'], it[1]]}.collect()
    // }
    // else if(bwa && "gridss" in callers) {
    //     ch_bwa_index_input = Channel.fromPath(bwa).map{[[id:"bwa"],it]}.collect()
    //     if(bwa.endsWith(".tar.gz")) {
    //         UNTAR_BWA(
    //             ch_bwa_index_input
    //         )
    //         ch_versions = ch_versions.mix(UNTAR_BWA.out.versions)

    //         UNTAR_BWA.out.untar
    //             .collect()
    //             .set { ch_bwa_index }
    //     } else {
    //         ch_bwa_index = ch_bwa_index_input
    //     }
    // }
    // else {
    //     ch_bwa_index = Channel.empty()
    // }

    if(annotate && !annotsv_annotations && callers.intersect(annotationCallers)) {
        ANNOTSV_INSTALLANNOTATIONS()
        ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)

        ANNOTSV_INSTALLANNOTATIONS.out.annotations
            .collect { annotations -> [[id:"annotsv_annotations"], annotations] }
            .set { ch_annotsv_annotations }
    }
    else if(annotate && callers.intersect(annotationCallers)) {
        ch_annotsv_annotations_input = Channel.fromPath(annotsv_annotations).collect { annotations -> [[id:"annotsv_annotations"], annotations] }
        if(annotsv_annotations.endsWith(".tar.gz")){
            UNTAR_ANNOTSV(
                ch_annotsv_annotations_input
            )
            ch_versions = ch_versions.mix(UNTAR_ANNOTSV.out.versions)

            UNTAR_ANNOTSV.out.untar
                .collect()
                .set { ch_annotsv_annotations }
        } else {
            ch_annotsv_annotations = Channel.fromPath(annotsv_annotations).collect { annotations -> [[id:"annotsv_annotations"], annotations] }
        }
    }
    else {
        ch_annotsv_annotations = Channel.empty()
    }

    if(!vep_cache && annotate && callers.intersect(annotationCallers)) {
        ENSEMBLVEP_DOWNLOAD(
            Channel.of([[id:"vep_cache"], vep_assembly, species, vep_cache_version]).collect()
        )
        ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

        ch_vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect { annotations -> annotations[1] }
    }
    else if (vep_cache && annotate && callers.intersect(annotationCallers)) {
        ch_vep_cache = Channel.fromPath(vep_cache).collect()
    }
    else {
        ch_vep_cache = Channel.empty()
    }

    //
    // Prepare the inputs
    //

    BAM_PREPARE_SAMTOOLS(
        ch_samplesheet.map { meta, cram, crai, small_variants ->
            [ meta, cram, crai ]
        },
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(BAM_PREPARE_SAMTOOLS.out.versions)

    ch_samplesheet
        .map { meta, cram, crai, small_variants ->
            [ meta, small_variants ]
        }
        .groupTuple()
        .map { meta, small_variants ->
            [ meta, small_variants.find { small_variant -> small_variant != [] } ?: [] ]
        }
        .set { ch_deduplicated }

    BAM_PREPARE_SAMTOOLS.out.crams
        .join(ch_deduplicated, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_input_no_sex }

    //
    // Determine the gender if needed
    //

    ch_input_sex = Channel.empty()
    if(callers.intersect(sexCallers)) {
        ch_input_no_sex
            .branch { meta, cram, crai, small_variants ->
                sex: meta.sex
                no_sex: !meta.sex
            }
            .set { ch_samplegender_input }

        NGSBITS_SAMPLEGENDER(
            ch_samplegender_input.no_sex,
            ch_fasta,
            ch_fai,
            "xy"
        )
        ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions.first())

        NGSBITS_SAMPLEGENDER.out.tsv
            .join(ch_samplegender_input.no_sex, failOnDuplicate:true, failOnMismatch:true)
            .map { meta, tsv, cram, crai ->
                def new_meta = meta + [sex:get_sex(tsv, meta.sample)]
                return [ new_meta, cram, crai ]
            }
            .mix(ch_samplegender_input.sex)
            .set { ch_input_sex }
    } else {
        ch_input_no_sex.set { ch_input_sex }
    }

    ch_input_sex
        .multiMap { meta, cram, crai, small_variants ->
            crams:          [ meta, cram, crai ]
            small_variants: [ meta, small_variants ]
        }
        .set { ch_inputs }

    //
    // Call the variants
    //

    if(sv_callers_to_use){

        count_types += 1
        variant_types.add("sv")

        BAM_SV_CALLING(
            ch_inputs.crams,
            ch_fasta,
            ch_fai,
            ch_manta_config,
            ch_svync_configs,
            sv_callers_to_use
        )

        ch_versions = ch_versions.mix(BAM_SV_CALLING.out.versions)
        ch_reports  = ch_reports.mix(BAM_SV_CALLING.out.reports)
        ch_annotation_input = ch_annotation_input.mix(BAM_SV_CALLING.out.vcfs)

    }

    //
    // Copy number calling
    //

    if(cnv_callers_to_use){

        count_types += 1
        variant_types.add("cnv")

        BAM_CNV_CALLING(
            ch_inputs.crams,
            ch_fasta,
            ch_fai,
            ch_qdnaseq_male,
            ch_qdnaseq_female,
            ch_wisecondorx_reference,
            ch_blacklist,
            ch_bedgovcf_configs,
            cnv_callers_to_use
        )
        ch_versions         = ch_versions.mix(BAM_CNV_CALLING.out.versions)
        ch_annotation_input = ch_annotation_input.mix(BAM_CNV_CALLING.out.vcfs)
    }

    //
    // Annotate the variants
    //

    if(annotate) {
        VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO(
            ch_annotation_input,
            ch_inputs.small_variants,
            ch_fasta,
            ch_fai,
            ch_annotsv_annotations,
            ch_annotsv_candidate_genes,
            ch_annotsv_gene_transcripts,
            ch_vep_cache,
            ch_vep_extra_files,
            ch_vcfanno_lua,
            val_vcfanno_resources,
            variant_types,
            genome,
            species,
            vep_cache_version,
            annotations_filter,
            val_vcfanno_toml,
            val_default_vcfanno_tomls
        )

        ch_reports  = ch_reports.mix(VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO.out.reports)
        ch_versions = ch_versions.mix(VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO.out.versions)
        ch_outputs  = ch_outputs.mix(VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO.out.vcfs)
    } else {
        ch_outputs  = ch_outputs.mix(ch_annotation_input)
    }

    //
    // Estimate repeat sizes
    //

    if(callers.intersect(repeatsCallers)){

        count_types += 1

        BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER(
            ch_inputs.crams,
            ch_fasta,
            ch_fai,
            ch_catalog
        )
        ch_versions = ch_versions.mix(BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER.out.versions)
        ch_outputs  = ch_outputs.mix(BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER.out.vcfs)

    }

    //
    // Concatenate the VCF files from different types of analysis
    //

    ch_concat_vcfs = Channel.empty()
    if(count_types > 1 && concat_output) {
        ch_outputs
            .map { meta, vcf, tbi ->
                def new_meta = meta - meta.subMap("variant_type")
                [ new_meta, vcf, tbi ]
            }
            .set { ch_concat_input }

        VCF_CONCAT_BCFTOOLS(
            ch_concat_input,
            count_types
        )
        ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

        VCF_CONCAT_BCFTOOLS.out.vcfs
            .set { ch_concat_vcfs }
    } else {
        ch_outputs.set { ch_concat_vcfs }
    }

    //
    // Merge VCFs of the same family into a multi-sample VCF
    //

    VCF_MERGE_FAMILY_JASMINE(
        ch_concat_vcfs,
        ch_fasta,
        ch_fai
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def get_sex(tsv, id) {
    if(workflow.stubRun) {
        return "female"
        log.warn("Couldn't define the sex of sample ${id}. Defaulting to female. (Specify the sex in the samplesheet to avoid this warning.)")
    }
    def split_tsv = tsv.splitCsv(sep:"\t", header:true, strip:true)
    def sex = split_tsv[0].gender
    if(sex == "others") {
        sex = "female"
        log.warn("Couldn't define the sex of sample ${id}. Defaulting to female. (Specify the sex in the samplesheet to avoid this warning.)")
    }
    return sex
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
