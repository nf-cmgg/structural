# nf-cmgg/structural: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.3.0dev

### `Added`

1. Added the samplesheet to the pipeline output as `OUTDIR/samplesheet.csv`
2. Added the `--bedpe` parameter. This makes the pipeline output BEDPE files alongside the VCF files.
3. Added parallelization on SV type to the delly flow
4. Added a `--gtf` parameter for annotation of gene and transcript overlap using `gatk SVAnnotate`.
5. Added `StrVCTVRE` as a new annotation tool

### `Changes`

1. Updated the pipeline to be compatible with future nextflow updates
2. Removed the `max_cpus`, `max_memory` and `max_time` parameter in favor of the new built-in `resourceLimits`
3. Replaced `nf-validation` with `nf-schema`
4. Updated to nf-core template v3.2.0
5. Fixed language server errors
6. Removed the old output publishing code and used the new workflow output definitions instead
7. Bumped the minimal nextflow version to 25.04.0
8. VCFanno will now run when `--vcfanno_toml` has been given and `--annotate` has not been given. You still need to supply `--annotate` to get the full annotation, but this can be used check for common variants without having to perform full annotation.
9. Changed the `--annotations_filter` parameter to a `--filter` parameter. This parameter takes an argument of `bcftools filter` to filter the resulting VCFs.
10. Removed the `--delly_sv_types` parameter.
11. Moved all `wisecondorx` and `qdnaseq` outputs to a separate directory in each sample output.
12. Bumped all annotation modules to the latest versions
13. Reworked the annotation structure to a per tool structure. Specify the annotations tools you want to run with `--annotate_tools`. This parameter takes a comma-separated list of tool names (options: `vep`, `vcfanno`, `svannotate`, `strvctvre` or `all` (=> all tools))

### `Fixed`

1. Fail the pipeline when the sex determination failed. This will now prompt the user to add the sex to the samplesheet so the pipeline doesn't do any wrong assumptions
2. Fixed the Jasmine module output VCFs being empty when no variants have been merged. This file now contains the header of one of the input VCFs
3. AnnotSV VCF files are now sorted before trying to combine it with the VEP output.
4. Fixed a map issue when the sex field is empty in the samplesheet

## v0.2.0 - [19 July 2024] - Mighty Manneken Pis

### `Added`

1. Added keyword shortcuts to the `--callers` parameter (these can also be used as comma-separated entries to the parameter):
   - `all`: Use all callers available in the pipeline
   - `sv`: Use all SV callers available in the pipeline
   - `cnv`: Use all CNV callers available in the pipeline
   - `rre`: Use all repeat region expansion callers available in the pipeline
2. Added the possibilty to annotate on HPO terms. Add the right HPO terms to the samplesheet in the `hpo` field

### `Changes`

1. Updated all WisecondorX modules to version 1.2.6 and added the `--seed` argument to `WisecondorX predict`
2. Removed support for the `phenotypes` VEP plugin. Common VEP plugin support will be added later
3. Made the main workflow pluggable, making it possible to use this pipeline in a meta pipeline
4. Updated all modules to their latest version

### `Fixed`

1. The smoove outputs are now correct when using `--output_callers`

## v0.1.0 - [3 April 2024] - Amazing Atomium

Initial release of nf-cmgg/structural, created with the [nf-core](https://nf-co.re/) template.
Initial release of nf-cmgg/structural, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
