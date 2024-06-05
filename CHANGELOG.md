# nf-cmgg/structural: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0dev

### `Added`

1. Added keyword shortcuts to the `--callers` parameter (these can also be used as comma-separated entries to the parameter):
   - `all`: Use all callers available in the pipeline
   - `sv`: Use all SV callers available in the pipeline
   - `cnv`: Use all CNV callers available in the pipeline
   - `rre`: Use all repeat region expansion callers available in the pipeline
2. Added the possibilty to annotate on HPO terms. Add the right HPO terms to the samplesheet in the `hpo` field

### `Changes`

1. Updated all WisecondorX modules to version 1.2.6 and added the `--seed` argument to `WisecondorX predict`
2. Removed support for the `phenotypes` VEP plugin. Commen VEP plugin support will be added later
3. Made the main workflow pluggable, making it possible to use this pipeline in a meta pipeline
4. Updated all modules to their latest version

### `Fixed`

1. The smoove outputs are now correct when using `--output_callers`

## v0.1.0 - [3 April 2024] - Amazing Atomium

Initial release of nf-cmgg/structural, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
