#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/structural
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/structural
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STRUCTURAL              } from './workflows/structural'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_structural_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_structural_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_structural_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Path to comma-separated file containing information about the samples in the experiment.
    input: String

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String

    // Email address for completion summary.
    email: String

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String

    // Name of iGenomes reference.
    genome: String

    // Path to FASTA genome file.
    fasta: String = getGenomeAttribute('fasta')

    // The index of the FASTA reference file
    fai: String = getGenomeAttribute('fai')

    // The sequence dictionary of the FASTA reference file
    dict: String = getGenomeAttribute('dict')

    // Path to GTF file for the reference genome. Gene and transcript annotations will be added when this file is provided
    gtf: String = getGenomeAttribute('gtf')

    // Path to the expansionhunter catalog
    expansionhunter_catalog: String = getGenomeAttribute('expansionhunter_catalog')

    // Path to the male qdnaseq reference file
    qdnaseq_male: String = getGenomeAttribute("qdnaseq_male_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString())

    // Path to the female qdnaseq reference file
    qdnaseq_female: String = getGenomeAttribute("qdnaseq_female_${params.qdnaseq_bin_size.toInteger() / 1000}kbp".toString())

    // Path to the wisecondorx reference file
    wisecondorx_reference: String? = getGenomeAttribute('wisecondorx_reference')

    // Path to the StrVCTVRE phylo bigwig file
    strvctvre_phylop: String? = getGenomeAttribute('strvctvre_phylop')

    // Path to the StrVCTVRE data directory
    strvctvre_data: String? = getGenomeAttribute('strvctvre_data')

    // Path to the blacklist BED file
    blacklist: String

    // Do not load the iGenomes reference config.
    igenomes_ignore: Boolean

    // The base path where the iGenomes references can be found
    igenomes_base: String

    // The base path where the references can be found
    genomes_base: String

    // Whether or not to use the references found in the `--genomes_base` folder
    genomes_ignore: Boolean

    // The config base path for the cmgg configs
    cmgg_config_base: String = '/conf/'

    // A map containing all references for all genomes
    genomes

    // Git commit id for Institutional configs.
    custom_config_version: String = 'master'

    // Base directory for Institutional configs.
    custom_config_base: String = 'https://raw.githubusercontent.com/nf-core/configs/master'

    // Institutional config name.
    config_profile_name: String

    // Institutional config description.
    config_profile_description: String

    // Institutional config contact information.
    config_profile_contact: String

    // Institutional config URL link.
    config_profile_url: String

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: String = 'https://raw.githubusercontent.com/nf-core/test-datasets/'

    // Display the help message.
    help

    // Display the full detailed help message.
    help_full: Boolean

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean

    // Display version and exit.
    version: Boolean

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String = 'copy'

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Incoming hook URL for messaging service
    hook_url: String

    // Custom config file to supply to MultiQC.
    multiqc_config: String

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: String

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: String

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Output monochrome logs
    monochrome_logs: Boolean

    trace_report_suffix: String

    // A comma-seperated list of callers to use. Can be one or more these: smoove/delly/manta/expansionhunter/qdnaseq/wisecondorx.
    callers: String = 'manta,smoove,delly,expansionhunter,wisecondorx'

    // Output the VCF files from different callers. Warning: This produces a lot of additional output and should only be used for testing purposes
    output_callers: Boolean

    // The minimum amount of SV callers that should detect a variant. All variants that have a lower amount of callers supporting it, will be removed. (Only used when more than one caller is given)
    sv_callers_support: Integer = 1

    // The minimum amount of CNV callers that should detect a variant. All variants that have a lower amount of callers supporting it, will be removed. (Only used when more than one caller is given)
    cnv_callers_support: Integer = 1

    // Run the annotation with Ensembl VEP and AnnotSV (and optionally VCFanno).
    annotate: Boolean

    // A comma-separated list of tools to use for annotation. Possible values: vep,svannotate,vcfanno,strvctvre. Default is all tools.
    annotate_tools: String = 'all'

    // Also output a concatenated VCF with all variant types analysed included.
    concat_output: Boolean

    // The filter options to perform on SV and CNV VCF files as postprocessing
    filter: String

    // Output BEDPE files derived from the VCF files alongside the VCF files
    bedpe: Boolean

    // The mapping quality to use for delly
    delly_map_qual: Integer = 1

    // The minimum clique size to use for delly
    delly_min_clique_size: Integer = 2

    // A config file to supply to manta
    manta_config: String

    // The bin size to use for qdnaseq.
    qdnaseq_bin_size: Integer = 100000

    // The minimum value of the absolute cnv ratio for a variant to be considered a CNV.
    qdnaseq_cnv_ratio: Float = 0.5

    // The genome assembly to download the cache of.
    vep_assembly: String = 'GRCh38'

    // The version of the VEP cache to use.
    vep_cache_version: Integer = 112

    // The path to the VEP cache folder
    vep_cache: String = getGenomeAttribute('vep_cache')

    // The version of VEP to use
    vep_version: Float = 112.0

    // The species used for the analysis. Should be all lowercase and spaces should be underscorses.
    species: String = 'homo_sapiens'

    // The full path to the VCFanno config TOML file. This file will be used to dynamically overwrite default configs for this pipeline run
    vcfanno_toml: String

    // The full path to a lua script for VCFanno
    vcfanno_lua: String

    // A comma-delimited list of files referenced in the VCFanno config and their indices.
    vcfanno_resources: String
}

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    def ch_multiqc_config          = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config   = params.multiqc_config ? channel.fromPath( params.multiqc_config, checkIfExists: true ) : channel.empty()
    def ch_multiqc_logo            = params.multiqc_logo   ? channel.fromPath( params.multiqc_logo, checkIfExists: true ) : channel.empty()
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    //
    // WORKFLOW: Run main workflow
    //
    STRUCTURAL (
        // channels
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_multiqc_config,
        ch_multiqc_custom_config,
        ch_multiqc_logo,
        ch_multiqc_custom_methods_description,

        // files
        params.fasta,
        params.fai,
        params.dict,
        params.gtf,
        params.expansionhunter_catalog ?: "https://github.com/Illumina/ExpansionHunter/raw/master/variant_catalog/grch38/variant_catalog.json",
        params.qdnaseq_female,
        params.qdnaseq_male,
        params.wisecondorx_reference,
        params.vep_cache,
        params.vcfanno_lua,
        params.vcfanno_resources,
        params.vcfanno_toml,
        params.blacklist,
        params.manta_config ?: "${projectDir}/assets/manta_config.ini",
        "${projectDir}/assets/svync",
        "${projectDir}/assets/bedgovcf",
        "${projectDir}/assets/vcfanno",
        params.strvctvre_phylop,
        params.strvctvre_data,

        // booleans
        params.annotate,
        params.concat_output,
        params.bedpe,

        // values
        params.callers,
        params.sv_callers_support,
        params.cnv_callers_support,
        params.genome,
        params.species,
        params.vep_assembly,
        params.vep_cache_version,
        params.filter,
        params.outdir,
        params.annotate_tools ? params.annotate_tools.tokenize(",") : []
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        STRUCTURAL.out.multiqc_report
    )

    publish:
    caller_vcfs     = STRUCTURAL.out.caller_vcfs
    sample_vcfs     = STRUCTURAL.out.sample_vcfs
    family_vcfs     = STRUCTURAL.out.family_vcfs
    qdnaseq_out     = STRUCTURAL.out.qdnaseq_out
    wisecondorx_out = STRUCTURAL.out.wisecondorx_out
    bedpe           = STRUCTURAL.out.bedpe
    multiqc         = STRUCTURAL.out.multiqc_report
    multiqc_data    = STRUCTURAL.out.multiqc_data

}

output {
    caller_vcfs {
        enabled params.output_callers
        path { meta, vcf, tbi ->
            vcf >> "${meta.sample}/${meta.caller}/${meta.sample}.vcf.gz"
            tbi >> "${meta.sample}/${meta.caller}/${meta.sample}.vcf.gz.tbi"
        }
    }
    sample_vcfs {
        path { meta, vcf, tbi ->
            vcf >> "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}" + ".vcf.gz"
            tbi >> "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}" + ".vcf.gz.tbi"
        }
    }
    family_vcfs {
        path { meta, vcf, tbi ->
            vcf >> "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}" + ".vcf.gz"
            tbi >> "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}" + ".vcf.gz.tbi"
        }
    }
    qdnaseq_out {
        path { meta, _bed_qdnaseq -> "$meta.id/qdnaseq/"
            // def base_suffix = bed_qdnaseq.name.replace(meta.id, "${meta.id}.qdnaseq")
            // bed_qdnaseq >> bed_qdnaseq.name == "statistics.out" ?
            //     "${meta.id}/${meta.id}.qdnaseq.statistics.out" :
            //     "${meta.id}/${base_suffix}"
        }
    }
    wisecondorx_out {
        path { meta, _bed -> "$meta.id/wisecondorx/"
            // if(bed.name.endsWith(".png")) {
            //     bed >> "${meta.id}/${meta.id}.wisecondorx.${bed.name}"
            // } else {
            //     def new_name = bed.name.replaceFirst(meta.id, "${meta.id}.wisecondorx")
            //     bed >> "${meta.id}/${new_name}"
            // }
        }
    }
    bedpe {
        path { meta, bedpe ->
            bedpe >> "${meta.id}/${meta.id}${meta.variant_type ? '.' + meta.variant_type : ''}" + ".bedpe"
        }
    }
    multiqc {
        path { "multiqc/" }
    }
    multiqc_data {
        path { "multiqc/" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
