# nf-cmgg/structural: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.3.0dev

### `Added`

### `Changes`

1. Updated the pipeline to be compatible with future nextflow updates
2. Removed the `max_cpus`, `max_memory` and `max_time` parameter in favor of the new built-in `resourceLimits`
3. Replaced `nf-validation` with `nf-schema`
4. Updated to nf-core template v3.2.0
5. Fixed language server errors

### `Fixed`

1. Fail the pipeline when the sex determination failed. This will now prompt the user to add the sex to the samplesheet so the pipeline doesn't do any wrong assumptions

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
