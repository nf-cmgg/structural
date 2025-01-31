# nf-cmgg/structural: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The files and directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [SV calling](#sv-calling) - Merged set of structural variant (SV) calls for the selected SV callers
- [CNV calling](#cnv-calling) - Merged set of copy number variant (CNV) calls for the selected CNV callers
- [RRE calling](#rre-calling) - Merged set of repeat region expansions (RRE) calls for the selected RRE callers
- [SV annotation](#sv-annotation) - Annotated set of structural variant (SV) calls for the selected SV callers
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### SV calling

<details markdown="1">
<summary>Output files</summary>

- `sampleID/sampleID.sv.vcf.gz`: vcf format file with merged SV calls for all selected callers.
- `sampleID/sampleID.sv.vcf.gz.tbi`: tabix index for the vcf format file with merged SV calls for all selected callers.

</details>

SV calling runs all selected callers individually and merges the calls after. Optionally, the output of the selected callers can also be stored individually. SV calling provides the following callers:

#### [Delly](https://github.com/dellytools/delly)

<details markdown="1">
<summary>Output files</summary>

- `sampleID/delly/`
  - `sampleID.delly.vcf.gz`: vcf format file with SV calls for Delly.
  - `sampleID.delly.vcf.gz.tbi`: tabix index for the vcf format file with SV calls for Delly.

</details>

[Delly](https://github.com/dellytools/delly) is an integrated structural variant (SV) prediction method that can discover, genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read and long-read massively parallel sequencing data. It uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome.

#### [Manta](https://github.com/Illumina/manta)

<details markdown="1">
<summary>Output files</summary>

- `sampleID/manta/`
  - `sampleID.manta.vcf.gz`: vcf format file with SV calls for Manta.
  - `sampleID.manta.vcf.gz.tbi`: tabix index for the vcf format file with SV calls for Manta.

 </details>

[Manta](https://github.com/Illumina/manta) calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs. Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a single efficient workflow. The method is designed for rapid analysis on standard compute hardware: NA12878 at 50x genomic coverage is analyzed in less than 20 minutes on a 20 core server, and most WGS tumor/normal analyses can be completed within 2 hours. Manta combines paired and split-read evidence during SV discovery and scoring to improve accuracy, but does not require split-reads or successful breakpoint assemblies to report a variant in cases where there is strong evidence otherwise. It provides scoring models for germline variants in small sets of diploid samples and somatic variants in matched tumor/normal sample pairs. There is experimental support for analysis of unmatched tumor samples as well. Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences in VCF 4.1 format. See the user guide for a full description of capabilities and limitations.

#### [Smoove](https://github.com/brentp/smoove)

<details markdown="1">
<summary>Output files</summary>

- `sampleID/smoove/`
  - `sampleID.smoove.vcf.gz`: vcf format file with SV calls for Smoove.
  - `sampleID.smoove.vcf.gz.tbi`: tabix index for the vcf format file with SV calls for Smoove.

 </details>

[Smoove](https://github.com/brentp/smoove) simplifies and speeds calling and genotyping SVs for short reads. It also improves specificity by removing many spurious alignment signals that are indicative of low-level noise and often contribute to spurious calls. It wraps existing software and adds some internal read-filtering to simplify calling and genotyping structural variants. It parallelizes each step as it can, for example, it streams lumpy output directly to multiple svtyper processes for genotyping.

### CNV calling

<details markdown="1">
<summary>Output files</summary>

- `sampleID/sampleID.cnv.vcf.gz`: vcf format file with merged CNV calls for all selected callers.
- `sampleID/sampleID.cnv.vcf.gz.tbi`: tabix index for the vcf format file with merged CNV calls for all selected callers.

</details>

CNV calling runs all selected callers individually and merges the calls after. Optionally, the output of the selected callers can also be stored individually. If this option is not selected, than files that would normally appear in the callers' result directories are provided in the main results directory. CNV calling provides the following callers:

#### [WisecondorX](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX)

<details markdown="1">
<summary>Output files</summary>

- `sampleID/wisecondorx/`
  - `sampleID.wisecondorx.vcf.gz`: vcf format file with CNV calls for WisecondorX.
  - `sampleID.wisecondorx.vcf.gz.tbi`: tabix index for the vcf format file with CNV calls for WisecondorX.
  - `sampleID.wisecondorx_aberrations.bed`: bed format file with aberrant segments.
  - `sampleID.wisecondorx_bins.bed`: bed format file with bin-wise information.
  - `sampleID.wisecondorx_segments.bed`: bed format file with segment-wise information.
  - `sampleID/chr1-X.png`: copy number profiles for every chromosome.
  - `sampleID/genome_wide.png`: genome-wide copy number profiles.

</details>

After extensively comparing different (shallow) whole-genome sequencing-based copy number detection tools, including WISECONDOR, QDNAseq, CNVkit, Control-FREEC, BIC-seq2 and cn.MOPS, WISECONDOR appeared to normalize sequencing data in the most consistent way, as shown by our paper. Nevertheless, WISECONDOR has limitations: Stouffer's Z-score approach is error-prone when dealing with large amounts of aberrations, the algorithm is extremely slow (24h) when operating at small bin sizes (15 kb), and sex chromosomes are not part of the analysis. Here, we present WisecondorX, an evolved WISECONDOR that aims at dealing with previous difficulties, resulting in overall superior results and significantly lower computing times, allowing daily diagnostic use. WisecondorX is meant to be applicable not only to NIPT, but also gDNA, PGT, FFPE, LQB, ... etc.

#### [QDNAseq](https://github.com/ccagc/QDNAseq)

<details markdown="1">
<summary>Output files</summary>

- `sampleID/qdnaseq/`
  - `sampleID.qdnaseq.abberations.bed`: bed format file with aberrant copy numbers.
  - `sampleID.qdnaseq.bed`: bed format file with copy numbers.
  - `sampleID.qdnaseq.cna`: file with bin-wise information.
  - `sampleID.qdnaseq.vcf.gz`: vcf format file with CNV calls for QDNAseq.
  - `sampleID.qdnaseq.vcf.gz.tbi`: tabix index for the vcf format file with CNV calls for QDNAseq.
  - `sampleID.qdnaseq_segments.txt`: file with segment-wise information.
  - `statistics.out`: statistics report.

</details>

Quantitative DNA sequencing for chromosomal aberrations. The genome is divided into non-overlapping fixed-sized bins, number of sequence reads in each counted, adjusted with a simultaneous two-dimensional loess correction for sequence mappability and GC content, and filtered to remove spurious regions in the genome. Downstream steps of segmentation and calling are also implemented via packages DNAcopy and CGHcall, respectively.

### RRE calling

<details markdown="1">
<summary>Output files</summary>

- `sampleID/sampleID.repeats.vcf.gz`: vcf format file with merged RRE calls for all selected callers.
- `sampleID/sampleID.repeats.vcf.gz.tbi`: tabix index for the vcf format file with merged RRE calls for all selected callers.

</details>

RRE calling runs all selected callers individually and merges the calls after. Optionally, the output of the selected callers can also be stored individually. CNV calling provides the following callers:

</details>

#### [ExpansionHunter](https://github.com/Illumina/ExpansionHunter)

<details markdown="1">
<summary>Output files</summary>

- `sampleID/expansionhunter/`
  - `sampleID/sampleID.expansionhunter.vcf.gz`: vcf format file with RRE calls for ExpansionHunter.
  - `sampleID/sampleID.expansionhunter.vcf.gz.tbi`: tabix index for the vcf format file with RRE calls for ExpansionHunter.

</details>

There are a number of regions in the human genome consisting of repetitions of short unit sequence (commonly a trimer). Such repeat regions can expand to a size much larger than the read length and thereby cause a disease. Fragile X Syndrome, ALS, and Huntington's Disease are well known examples. [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) aims to estimate sizes of such repeats by performing a targeted search through a BAM/CRAM file for reads that span, flank, and are fully contained in each repeat.

### SV annotation

<details markdown="1">
<summary>Output files</summary>

- `sampleID/sampleID.sv.annotated.vcf.gz`: vcf format file with merged and annotated SV calls for all selected callers.
- `sampleID/sampleID.sv.annotated.vcf.gz.tbi`: tabix index for the vcf format file with merged and annotated SV calls for all selected callers.

</details>

SV annotation runs each annotation method individually and merges their output after. SV calling provides the following annotation methods:

#### [AnnotSV](https://github.com/lgmgeo/AnnotSV)

</details>

[AnnotSV](https://github.com/lgmgeo/AnnotSV) is a program designed for annotating and ranking Structural Variations (SV). This tool compiles functionally, regulatory and clinically relevant information and aims at providing annotations useful to i. interpret SV potential pathogenicity and ii. filter out SV potential false positives. Different types of SV exist including deletions, duplications, insertions, inversions, translocations or more complex rearrangements. They can be either balanced or unbalanced. When unbalanced and resulting in a gain or loss of material, they are called Copy Number Variations (CNV). CNV can be described by coordinates on one chromosome, with the start and end positions of the SV (deletions, insertions, duplications).Complex rearrangements with several breakends can arbitrary be summarized as a set of novel adjacencies, as described
in the Variant Call Format specification VCFv4.3.

#### [ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)

</details>

[ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) determines the effect of your variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
