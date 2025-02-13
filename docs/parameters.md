# nf-cmgg/structural pipeline parameters

A bioinformatics best-practice analysis pipeline for calling structural variants (SVs), copy number variants (CNVs) and repeat region expansions (RREs) from short DNA reads

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>| `string` |  | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `genome` | Name of iGenomes reference. <details><summary>Help</summary><small>If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`. <br><br>See the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details.</small></details>| `string` |  |  |  |
| `fasta` | Path to FASTA genome file. <details><summary>Help</summary><small>This parameter is *mandatory* if `--genome` is not specified. </small></details>| `string` |  |  |  |
| `fai` | The index of the FASTA reference file | `string` |  |  |  |
| `expansionhunter_catalog` | Path to the expansionhunter catalog | `string` |  |  |  |
| `qdnaseq_male` | Path to the male qdnaseq reference file | `string` |  |  |  |
| `qdnaseq_female` | Path to the female qdnaseq reference file | `string` |  |  |  |
| `wisecondorx_reference` | Path to the wisecondorx reference file | `string` |  |  |  |
| `blacklist` | Path to the blacklist BED file | `string` |  |  |  |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>| `boolean` |  |  | True |
| `igenomes_base` | The base path where the iGenomes references can be found | `string` |  |  |  |
| `genomes_base` | The base path where the references can be found | `string` |  |  |  |
| `genomes_ignore` | Whether or not to use the references found in the `--genomes_base` folder | `boolean` |  |  |  |
| `cmgg_config_base` | The config base path for the cmgg configs | `string` | /conf/ |  |  |
| `genomes` | A map containing all references for all genomes | `object` |  |  | True |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display version and exit. | `boolean` |  |  | True |
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `monochrome_logs` | Output monochrome logs | `boolean` |  |  |  |
| `pipelines_testdata_base_path` | The base path to of the pipeline test data | `string` | https://raw.githubusercontent.com/nf-core/test-datasets/ |  | True |
| `trace_report_suffix` |  | `string` |  |  | True |

## Pipeline specific options

Options specific to the execution of this pipeline

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `callers` | A comma-seperated list of callers to use. Can be one or more these: smoove/delly/manta/expansionhunter/qdnaseq/wisecondorx. <details><summary>Help</summary><small>Following keywords can also be used as shortcuts to certain callers: <br>- all: Use all callers<br>- sv: Use all SV callers<br>- cnv: Use all CNV callers<br>- rre: Use all RRE callers</small></details>| `string` | manta,smoove,delly,expansionhunter,wisecondorx |  |  |
| `output_callers` | Output the VCF files from different callers. Warning: This produces a lot of additional output and should only be used for testing purposes | `boolean` |  |  |  |
| `sv_callers_support` | The minimum amount of SV callers that should detect a variant. All variants that have a lower amount of callers supporting it, will be removed. (Only used when more than one caller is given) | `integer` | 1 |  |  |
| `cnv_callers_support` | The minimum amount of CNV callers that should detect a variant. All variants that have a lower amount of callers supporting it, will be removed. (Only used when more than one caller is given) | `integer` | 1 |  |  |
| `annotate` | Run the annotation with Ensembl VEP and AnnotSV (and optionally VCFanno). | `boolean` |  |  |  |
| `concat_output` | Also output a concatenated VCF with all variant types analysed included. | `boolean` |  |  |  |
| `filter` | The filter options to perform on SV and CNV VCF files as postprocessing <details><summary>Help</summary><small>A common use case for this is to filter out common variants. You can supply a VCF file with common variants to VCFanno and filter out the matching variants using this parameter</small></details>| `string` |  |  |  |

## Delly parameters

Options specific for the Delly execution

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `delly_sv_types` | Which SV types delly should search for in the variant calling | `string` | ALL |  |  |
| `delly_map_qual` | The mapping quality to use for delly | `integer` | 1 |  |  |
| `delly_min_clique_size` | The minimum clique size to use for delly | `integer` | 2 |  |  |

## Manta parameters

Options specific for the Manta execution

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `manta_config` | A config file to supply to manta | `string` |  |  |  |

## qDNAseq parameters

Options specific for the qDNAseq execution

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `qdnaseq_bin_size` | The bin size to use for qdnaseq. | `integer` | 100000 |  |  |
| `qdnaseq_cnv_ratio` | The minimum value of the absolute cnv ratio for a variant to be considered a CNV. | `number` | 0.5 |  |  |

## VEP options

Options for the annotation with VEP

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `vep_assembly` | The genome assembly to download the cache of. <details><summary>Help</summary><small>This only needs to be supplied when no vep_cache is given</small></details>| `string` | GRCh38 |  |  |
| `vep_cache_version` | The version of the VEP cache to use. <details><summary>Help</summary><small>This version should be present in the folder supplied by `--vep_cache`. This version should be the same as `--vep_version` when no VEP cache was given with `--vep_cache`</small></details>| `integer` | 112 |  |  |
| `vep_cache` | The path to the VEP cache folder | `string` |  |  |  |
| `vep_version` | The version of VEP to use | `number` | 112.0 |  |  |
| `species` | The species used for the analysis. Should be all lowercase and spaces should be underscorses. | `string` | homo_sapiens |  |  |

## AnnotSV options

Options specific to the execution of AnnotSV

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `annotsv_annotations` | The full path to the AnnotSV annotations folder. This can be a tarzipped folder. <details><summary>Help</summary><small>When --annotate is set to true and this isn't given, the annotations will be downloaded automatically (this is not recommended though).</small></details>| `string` |  |  |  |
| `annotsv_candidate_genes` | The full path the candidate genes file for AnnotSV | `string` |  |  |  |
| `annotsv_gene_transcripts` | The full path to the gene transcripts file for AnnotSV | `string` |  |  |  |

## VCFanno options

Options for the execution of AnnotSV

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `vcfanno_toml` | The full path to the VCFanno config TOML file. This file will be used to dynamically overwrite default configs for this pipeline run | `string` |  |  |  |
| `vcfanno_lua` | The full path to a lua script for VCFanno | `string` |  |  |  |
| `vcfanno_resources` | A comma-delimited list of files referenced in the VCFanno config and their indices. <details><summary>Help</summary><small>This can contain glob patterns as well. Place all filenames together between double qoutes to not cause any irregularities</small></details>| `string` |  |  |  |
