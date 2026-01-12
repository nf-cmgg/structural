/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
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
include { BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER } from '../subworkflows/local/bam_repeat_estimation_expansionhunter/main'
include { VCF_ANNOTATE                          } from '../subworkflows/local/vcf_annotate/main'
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
include { GATK4_CREATESEQUENCEDICTIONARY    } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { PREPROCESS_GTF                    } from '../modules/local/preprocess_gtf/main'
include { BWA_INDEX                         } from '../modules/nf-core/bwa/index/main'
include { ENSEMBLVEP_DOWNLOAD               } from '../modules/nf-core/ensemblvep/download/main'
include { UNTAR as UNTAR_ANNOTSV            } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BWA                } from '../modules/nf-core/untar/main'
include { NGSBITS_SAMPLEGENDER              } from '../modules/nf-core/ngsbits/samplegender/main'
include { BCFTOOLS_FILTER                   } from '../modules/nf-core/bcftools/filter/main'
include { SVTOOLS_VCFTOBEDPE                } from '../modules/nf-core/svtools/vcftobedpe/main'
include { GATK4_SVANNOTATE                  } from '../modules/nf-core/gatk4/svannotate/main'
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
    dict                        // The dictionary of the fasta file
    gtf                         // The GTF file to use for annotation
    expansionhunter_catalog     // The expansionhunter catalog
    qdnaseq_female              // The QDNAseq annotations for female samples
    qdnaseq_male                // The QDNAseq annotations for male samples
    wisecondorx_reference       // The WisecondorX annotations file
    vep_cache                   // The VEP cache directory
    vcfanno_lua                 // A Lua script to use with vcfanno
    vcfanno_resources           // A comma delimited list of paths to vcfanno resource files
    vcfanno_toml                // The vcfanno config to
    blacklist                   // A BED file of blacklisted regions
    manta_config                // A configuration file to be used in Manta
    svync_dir                   // A directory containing svync configs (they must end in '.yaml' and contain the full name of the caller)
    bedgovcf_dir                // A directory containing bedgovcf configs (they must end in '.yaml' and contain the full name of the caller)
    vcfanno_default_dir         // A directory containing the default vcfanno configs (they must end in '.toml')
    strvctvre_phylop            // A bigwig file containing the phylop reference for StrVCTVRE
    strvctvre_data              // A directory containing the reference data for StrVCTVRE

    // boolean inputs
    annotate                    // Run annotation on SV and CNV data
    concat_output               // Concatenate all output files
    bedpe                       // Create BEDPE files from the VCFs and output them

    // value inputs
    input_callers               // All callers to be used
    sv_callers_support          // How many SV callers are needed to support the variant
    cnv_callers_support         // How many CNV callers are needed to support the variant
    genome                      // The genome to use
    species                     // The species to be used by VEP
    vep_assembly                // The genome assembly to be downloaded for VEP
    vep_cache_version           // The version of the VEP cache to use
    filter                      // The filter pattern to use after annotation
    outdir                      // The output directory of the pipeline
    annotate_tools              // The tools to be used for annotation

    main:

    def ch_versions         = channel.empty()
    def ch_reports          = channel.empty()
    def ch_caller_vcfs      = channel.empty()
    def ch_multiqc_files    = channel.empty()

    def variant_types = [] // The variant types that can be annotated this run
    def count_types = 0 // The amount of different variant types that can be concatenated

    //
    // Create input channels from parameters
    //

    def ch_fasta                    = channel.fromPath(fasta).collect { fasta_file -> [[id:'fasta'], fasta_file ] }
    // def ch_annotsv_candidate_genes  = annotsv_candidate_genes ?  channel.fromPath(annotsv_candidate_genes).collect { genes_file -> [[], genes_file] } : [[],[]]
    // def ch_annotsv_gene_transcripts = annotsv_gene_transcripts ? channel.fromPath(annotsv_gene_transcripts).collect { transcripts_file -> [[], transcripts_file] } : [[],[]]
    def ch_vcfanno_lua              = vcfanno_lua ?              channel.fromPath(vcfanno_lua).collect() : []
    def ch_catalog                  = expansionhunter_catalog ?  channel.fromPath(expansionhunter_catalog).collect { catalog_file -> [[id:'catalog'], catalog_file] } : [[],[]]
    def ch_qdnaseq_male             = qdnaseq_male ?             channel.fromPath(qdnaseq_male).collect { qdnaseq_file -> [[id:'qdnaseq_male'], qdnaseq_file] } : [[],[]]
    def ch_qdnaseq_female           = qdnaseq_female ?           channel.fromPath(qdnaseq_female).collect { qdnaseq_file -> [[id:'qdnaseq_female'], qdnaseq_file] } : [[],[]]
    def ch_wisecondorx_reference    = wisecondorx_reference ?    channel.fromPath(wisecondorx_reference).collect { wcx_file -> [[id:'wisecondorx'], wcx_file] } : [[],[]]
    def ch_blacklist                = blacklist ?                channel.fromPath(blacklist).collect { blacklist_file -> [[id:'blacklist'], blacklist_file] } : [[],[]]
    def ch_manta_config             = manta_config ?             channel.fromPath(manta_config).collect() : [[]]
    def ch_svync_configs            = svync_dir ?                channel.fromPath("${svync_dir}/*.yaml", checkIfExists:true).collect() : []
    def ch_bedgovcf_configs         = bedgovcf_dir ?             channel.fromPath("${bedgovcf_dir}/*.yaml", checkIfExists:true).collect() : []
    def ch_strvctvre_phylop         = strvctvre_phylop ?         channel.value([[id:'phylop'], file(strvctvre_phylop)]) : [[:],[]]
    def ch_strvctvre_data           = strvctvre_data ?           channel.value([[id:'strvctvre_data'], file(strvctvre_data)]) : [[:],[]]

    def val_vcfanno_resources       = vcfanno_resources ?        vcfanno_resources.split(",").collect { resource_file -> file(resource_file, checkIfExists:true) }.flatten() : []
    def val_default_vcfanno_tomls   = vcfanno_default_dir ?      files("${vcfanno_default_dir}/*.toml", checkIfExists:true) : []
    def val_vcfanno_toml            = vcfanno_toml ?             file(vcfanno_toml, checkIfExists:true) : []

    def ch_vep_extra_files = []

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
    def repeats_callers_to_use = callers.intersect(repeatsCallers)

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

    if (repeats_callers_to_use.size() > 0 && bedpe && concat_output) {
        error("Can't create BEDPE files from VCFs that contains repeat expansions. Don't specify either --concat_output or omit all repeat callers from the --callers parameter.")
    }

    if (annotate && (annotate_tools.contains("svannotate") || annotate_tools.contains("all")) && !gtf) {
        error("The GTF file is required when using SVAnnotate. Please provide it using the 'gtf' parameter.")
    }

    //
    // Create optional inputs
    //

    def ch_fai = channel.empty()
    if(!fai){
        SAMTOOLS_FAIDX(
            ch_fasta,
            [[], []]
        )

        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_fai      = SAMTOOLS_FAIDX.out.fai.collect { fai_file -> [[id:'fai'], fai_file] }
    }
    else {
        ch_fai = channel.fromPath(fai).collect { fai_file -> [[id:'fai'], fai_file] }
    }

    def ch_dict = channel.empty()
    // Dictionary is only needed for GATK4_SVANNOTATE
    if(!dict && gtf) {
        GATK4_CREATESEQUENCEDICTIONARY(
            ch_fasta
        )
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

        ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
    }
    else if(dict) {
        ch_dict = channel.fromPath(dict).collect { dict_file -> [[id:'dict'], dict_file] }
    }

    def ch_preprocessed_gtf = channel.empty()
    // Sanitize GTF file to adhere to the extremely strict GTF parsing in SVAnnotate
    if (gtf) {
        ch_sanitize_input = channel.fromPath(gtf).collect { gtf_file -> [[id:'gtf'], gtf_file] }
        PREPROCESS_GTF(
            ch_sanitize_input
        )
        ch_versions = ch_versions.mix(PREPROCESS_GTF.out.versions)

        ch_preprocessed_gtf = PREPROCESS_GTF.out.gtf.collect()
    }

    def ch_vep_cache = channel.empty()
    if(!vep_cache && annotate && callers.intersect(annotationCallers)) {
        ENSEMBLVEP_DOWNLOAD(
            channel.of([[id:"vep_cache"], vep_assembly, species, vep_cache_version]).collect()
        )
        ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

        ch_vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache.collect { annotations -> annotations[1] }
    }
    else if (vep_cache && annotate && callers.intersect(annotationCallers)) {
        ch_vep_cache = channel.fromPath(vep_cache).collect()
    }

    //
    // Prepare the inputs
    //

    BAM_PREPARE_SAMTOOLS(
        ch_samplesheet,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(BAM_PREPARE_SAMTOOLS.out.versions)

    def ch_input_no_sex = BAM_PREPARE_SAMTOOLS.out.crams

    //
    // Determine the gender if needed
    //

    def ch_inputs = channel.empty()
    if(callers.intersect(sexCallers)) {
        def ch_samplegender_input = ch_input_no_sex
            .branch { meta, _cram, _crai ->
                sex: meta.sex
                no_sex: !meta.sex
            }

        NGSBITS_SAMPLEGENDER(
            ch_samplegender_input.no_sex,
            ch_fasta,
            ch_fai,
            "xy"
        )
        ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions.first())

        ch_inputs = NGSBITS_SAMPLEGENDER.out.tsv
            .join(ch_samplegender_input.no_sex, failOnDuplicate:true, failOnMismatch:true)
            .map { meta, tsv, cram, crai ->
                def new_meta = meta + [sex:get_sex(tsv, meta.sample)]
                return [ new_meta, cram, crai ]
            }
            .mix(ch_samplegender_input.sex)
    } else {
        ch_inputs = ch_input_no_sex
    }

    //
    // Call the variants
    //

    def ch_annotation_input = channel.empty()
    if(sv_callers_to_use){

        count_types += 1
        variant_types.add("sv")

        BAM_SV_CALLING(
            ch_inputs,
            ch_fasta,
            ch_fai,
            ch_manta_config,
            ch_svync_configs,
            sv_callers_to_use
        )

        ch_caller_vcfs = ch_caller_vcfs.mix(BAM_SV_CALLING.out.caller_vcfs)
        ch_versions = ch_versions.mix(BAM_SV_CALLING.out.versions)
        ch_reports  = ch_reports.mix(BAM_SV_CALLING.out.reports)
        ch_annotation_input = ch_annotation_input.mix(BAM_SV_CALLING.out.vcfs)

    }

    //
    // Copy number calling
    //

    def ch_wisecondorx_out = channel.empty()
    def ch_qdnaseq_out = channel.empty()
    if(cnv_callers_to_use){

        count_types += 1
        variant_types.add("cnv")

        BAM_CNV_CALLING(
            ch_inputs,
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
        ch_wisecondorx_out  = BAM_CNV_CALLING.out.wisecondorx
        ch_qdnaseq_out      = BAM_CNV_CALLING.out.qdnaseq
    }

    //
    // Annotate the variants
    //

    def ch_annotation_output = ch_annotation_input
    if(annotate && annotate_tools.size() > 0) {
        VCF_ANNOTATE(
            ch_annotation_input,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_preprocessed_gtf,
            ch_vep_cache,
            ch_vep_extra_files,
            ch_vcfanno_lua,
            val_vcfanno_toml,
            ch_strvctvre_phylop,
            ch_strvctvre_data,
            genome,
            species,
            vep_cache_version,
            val_vcfanno_resources,
            val_default_vcfanno_tomls,
            annotate_tools
        )
        ch_versions = ch_versions.mix(VCF_ANNOTATE.out.versions)
        ch_reports  = ch_reports.mix(VCF_ANNOTATE.out.reports)
        ch_annotation_output = VCF_ANNOTATE.out.vcfs
    }

    def ch_outputs = ch_annotation_output
    if(filter) {
        BCFTOOLS_FILTER(
            ch_annotation_output
        )
        ch_outputs = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions)
    }

    //
    // Estimate repeat sizes
    //

    if(callers.intersect(repeatsCallers)){

        count_types += 1

        BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER(
            ch_inputs,
            ch_fasta,
            ch_fai,
            ch_catalog
        )

        ch_caller_vcfs  = ch_caller_vcfs.mix(BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER.out.caller_vcfs)
        ch_versions     = ch_versions.mix(BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER.out.versions)
        ch_outputs      = ch_outputs.mix(BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER.out.vcfs)

    }

    //
    // Concatenate the VCF files from different types of analysis
    //

    def ch_concat_vcfs = channel.empty()
    if(count_types > 1 && concat_output) {
        def ch_concat_input = ch_outputs
            .map { meta, vcf, tbi ->
                def new_meta = meta - meta.subMap("variant_type", "caller")
                [ new_meta, vcf, tbi ]
            }

        VCF_CONCAT_BCFTOOLS(
            ch_concat_input,
            count_types
        )
        ch_versions = ch_versions.mix(VCF_CONCAT_BCFTOOLS.out.versions)

        ch_concat_vcfs = VCF_CONCAT_BCFTOOLS.out.vcfs
    } else {
        ch_concat_vcfs = ch_outputs
    }

    //
    // Merge VCFs of the same family into a multi-sample VCF
    //

    def ch_merge_vcfs = ch_concat_vcfs.filter { meta, _vcf, _tbi ->
        meta.family_count > 1
    }

    VCF_MERGE_FAMILY_JASMINE(
        ch_merge_vcfs,
        ch_fasta,
        ch_fai
    )

    def ch_family_vcfs = VCF_MERGE_FAMILY_JASMINE.out.vcfs

    //
    // Convert VCF files to BEDPE files
    //

    def ch_bedpe = channel.empty()
    if(bedpe) {
        def ch_vcftobedpe_input = ch_concat_vcfs
            .mix(ch_family_vcfs)
            .filter { meta, _vcf, _tbi ->
                meta.variant_type != "repeats"
            }
            .map { meta, vcf, _tbi ->
                [ meta, vcf ]
            }

        SVTOOLS_VCFTOBEDPE(
            ch_vcftobedpe_input
        )
        ch_versions = ch_versions.mix(SVTOOLS_VCFTOBEDPE.out.versions)
        ch_bedpe = SVTOOLS_VCFTOBEDPE.out.bedpe
    }

    //
    // Collate and save software versions
    //

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'structural_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    caller_vcfs     = ch_caller_vcfs              // channel: [ val(meta), path(vcf), path(tbi) ]
    sample_vcfs     = ch_concat_vcfs              // channel: [ val(meta), path(vcf), path(tbi) ]
    family_vcfs     = ch_family_vcfs              // channel: [ val(meta), path(vcf), path(tbi) ]
    qdnaseq_out     = ch_qdnaseq_out              // channel: [ val(meta), path(file) ]
    wisecondorx_out = ch_wisecondorx_out          // channel: [ val(meta), path(file) ]
    bedpe           = ch_bedpe                    // channel: [ val(meta), path(bedpe) ]
    multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    multiqc_data    = MULTIQC.out.data            // channel: /path/to/multiqc_data
    versions        = ch_versions                 // channel: [ path(versions.yml) ]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def get_sex(tsv, id) {
    if(workflow.stubRun) {
        return "female"
        log.warn("STUB: Couldn't define the sex of sample ${id}. Defaulting to female. (Specify the sex in the samplesheet to avoid this warning.)")
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
