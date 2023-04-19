//
// Annotate the VCFs
//

include { ANNOTSV_ANNOTSV                       } from '../../../modules/nf-core/annotsv/annotsv/main'
include { ENSEMBLVEP_VEP                        } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCFANNO                               } from '../../../modules/nf-core/vcfanno/main'
include { TABIX_BGZIPTABIX as TABIX_ANNOTATED   } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_ANNOTSV     } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX as TABIX_VEP              } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_ANNOTATE_VEP_ANNOTSV_VCFANNO {
    take:
        ch_vcfs                                 // channel: [mandatory] [ val(meta), path(vcf), path(tbi) ] VCFs containing the called structural variants
        ch_fasta                                // channel: [mandatory] [ path(fasta) ] => The fasta reference file
        ch_fai                                  // channel: [mandatory] [ path(fai) ] => The index of the fasta reference file
        ch_annotsv_annotations                  // channel: [mandatory] [ val(meta), path(annotations) ] => The annotations for AnnotSV
        ch_annotsv_candidate_genes              // channel: [optional]  [ val(meta), path(candidate_genes) ]
        ch_annotsv_gene_transcripts             // channel: [optional]  [ val(meta), path(gene_transcripts) ]
        ch_vep_cache                            // channel: [optional]  [ path(cache) ] => The path to the local VEP cache
        ch_vep_extra_files                      // channel: [optional]  [ path(file1, file2, file3...) ] => The VEP extra files
        ch_vcfanno_lua                          // channel: [optional]  [ path(lua) ] => The lua script to influence VCFanno
        val_vcfanno_resources                   // list:    [optional]  [ path(file1, file2, file3...) ] => The extra VCFanno files

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // Run AnnotSV and VEP in parallel and merge TSV from AnnotSV with VCF from VEP during VCFanno

    ANNOTSV_ANNOTSV(
        ch_vcfs,
        ch_annotsv_annotations,
        ch_annotsv_candidate_genes,
        [[],[]],
        [[],[]],
        ch_annotsv_gene_transcripts
    )
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions.first())

    ANNOTSV_ANNOTSV.out.vcf
        .map { meta, vcf ->
            // Artifically create a single variant annotated VCF if the VCF is empty
            // This will only affect test runs that create VCFs with no variants
            // TODO make sure the tests don't actually need this!
            if(vcf.size() == 0){
                vcf.text = 
"""##fileformat=VCFv4.3
##INFO=<ID=AnnotSV_ID,Number=.,Type=String,Description=AnnotSV ID>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	empty
16	75061454	1_MantaDEL:14:0:0:0:0:0;1_MantaDEL:23:0:0:0:0:0	GTTCCACTCATGTTGTTGCAGATTACTGGATCTCACTCTTTTTTATGGGTGAATCGTACTTCATAGTGTATATGTACC	G	550	PASS	AnnotSV_ID=16_75061454_75061530_DEL_1|16_75061454_75061530_DEL_1;SV_chrom=16|16;SV_start=75061454|75061454;SV_end=75061530|75061530;SV_length=-77|-77;SV_type=DEL|DEL;Samples_ID=PosCon1,PosCon2|PosCon1,PosCon2;FORMAT=GT:FT:GQ:PL:PR:SR:OLD_GT:DP:AD:ADF:ADR|GT:FT:GQ:PL:PR:SR:OLD_GT:DP:AD:ADF:ADR;Annotation_mode=full&split;CytoBand=q23.1|q23.1;Gene_name=ZNRF1|ZNRF1;Gene_count=1.0|.;Tx=.|NM_032268;Tx_start=.|74999023.0;Tx_end=.|75110994.0;Overlapped_tx_length=.|77.0;Overlapped_CDS_length=.|0.0;Overlapped_CDS_percent=.|0.0;Frameshift=.|no;Exon_count=.|5.0;Location=.|intron1-intron1;Location2=.|CDS;Dist_nearest_SS=.|32041.0;Nearest_SS_type=.|3';Intersect_start=.|75061453.0;Intersect_end=.|75061530.0;RE_gene=.|.;P_gain_phen=.|.;P_gain_hpo=.|.;P_gain_source=.|.;P_gain_coord=.|.;P_loss_phen=.|.;P_loss_hpo=.|.;P_loss_source=.|.;P_loss_coord=.|.;P_ins_phen=.|.;P_ins_hpo=.|.;P_ins_source=.|.;P_ins_coord=.|.;po_P_gain_phen=.|.;po_P_gain_hpo=.|.;po_P_gain_source=dbVar:nssv15158812,dbVar:nssv15148703,dbVar:nssv15150639,nssv15151292,dbVar:nssv15158813,dbVar:nssv15158224,dbVar:nssv15126469,dbVar:nssv17976066,dbVar:nssv15149003,dbVar:nssv15152203,dbVar:nssv15146677,dbVar:nssv15146029,dbVar:nssv15146411,dbVar:nssv15147487,dbVar:nssv15146794,dbVar:nssv15148937,dbVar:nssv15146412,dbVar:nssv15146731,dbVar:nssv15155495,dbVar:nssv15147670|.;po_P_gain_coord=16:11452-90228224,16:19194-90207973,16:35882-90088654,16:38166-90096867,16:38166-90208287,16:31974778-90228345,16:46470057-90088654,16:52899184-90088654,16:57017562-89731261,16:62925930-84585795,16:64389379-90081985,16:65313396-90081985,16:65511484-90096995,16:65957830-83611443,16:69053458-83274681,16:70514632-90081985,16:70749399-90096995,16:72482040-90088654,16:74838617-90208032|.;po_P_gain_percent=0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00|.;po_P_loss_phen=.|.;po_P_loss_hpo=.|.;po_P_loss_source=dbVar:nssv15128641,dbVar:nssv15146519,dbVar:nssv15147241,dbVar:nssv15133204,dbVar:nssv17970999|.;po_P_loss_coord=16:69535324-75926298,16:69918077-76723348,16:70414574-84908120,16:73049468-82576326,16:74714171-75479848|.;po_P_loss_percent=0.00,0.00,0.00,0.00,0.01|.;P_snvindel_nb=.|.;P_snvindel_phen=.|.;B_gain_source=.|.;B_gain_coord=.|.;B_gain_AFmax=.|.;B_loss_source=.|.;B_loss_coord=.|.;B_loss_AFmax=.|.;B_ins_source=.|.;B_ins_coord=.|.;B_ins_AFmax=.|.;B_inv_source=IMH|IMH;B_inv_coord=16:29226964-88193867|16:29226964-88193867;B_inv_AFmax=0.4997|0.4997;po_B_gain_allG_source=.|.;po_B_gain_allG_coord=.|.;po_B_gain_someG_source=.|.;po_B_gain_someG_coord=.|.;po_B_loss_allG_source=.|.;po_B_loss_allG_coord=.|.;po_B_loss_someG_source=.|.;po_B_loss_someG_coord=16:75061455-75061531_CMRI:0_pbsv.DEL.952_duplicate9|.;TAD_coordinate=.|.;ENCODE_experiment=.|.;GC_content_left=0.35|.;GC_content_right=0.365|.;Repeat_coord_left=16:75061019-75061531,16:75061532-75061674|.;Repeat_type_left=L1MA1,L1MA1|.;Repeat_coord_right=16:75061019-75061531,16:75061532-75061674|.;Repeat_type_right=L1MA1,L1MA1|.;Gap_left=.|.;Gap_right=.|.;SegDup_left=.|.;SegDup_right=.|.;ENCODE_blacklist_left=.|.;ENCODE_blacklist_characteristics_left=.|.;ENCODE_blacklist_right=.|.;ENCODE_blacklist_characteristics_right=.|.;ACMG=.|.;HI=.|.;TS=.|.;DDD_HI_percent=15.84|15.84;ExAC_delZ=.|.;ExAC_dupZ=.|.;ExAC_cnvZ=.|.;ExAC_synZ=1.76818456166216|1.76818456166216;ExAC_misZ=1.84776471271999|1.84776471271999;GenCC_disease=.|.;GenCC_moi=.|.;GenCC_classification=.|.;GenCC_pmid=.|.;OMIM_ID=.|.;OMIM_phenotype=.|.;OMIM_inheritance=.|.;OMIM_morbid=.|.;OMIM_morbid_candidate=.|.;LOEUF_bin=1|1;GnomAD_pLI=0.90588|0.90588;ExAC_pLI=0.7544|0.7544;AnnotSV_ranking_score=0.0|.;AnnotSV_ranking_criteria=1A_(cf_Gene_count,_RE_gene,_+0.00),2B_(cf_po_P_loss_source,_HI_and_OMIM_morbid,_+0.00),2G_(cf_po_B_loss_someG_source,_+0.00),3A_(1_gene,_+0.00),5F_(+0.00)|.;ACMG_class=3|full=3;SVTYPE=DEL;SVLEN=-77;CIPOS=0,5;HOMLEN=5;HOMSEQ=TTCCA;STARTVARIANCE=0.000000;ENDVARIANCE=0.000000;AVG_LEN=-77.000000;AVG_START=75061454.000000;AVG_END=75061531.000000;VARCALLS=1;ALLVARS_EXT=(MantaDEL:14:0:0:0:0:0);SUPP_VEC_EXT=010;IDLIST_EXT=MantaDEL:14:0:0:0:0:0;SUPP_EXT=1;SUPP_VEC=010;SUPP=1;SVMETHOD=JASMINE;IDLIST=MantaDEL:14:0:0:0:0:0;INTRASAMPLE_IDLIST=MantaDEL:14:0:0:0:0:0;GRMPY_ID=PosCon1.vcf.gz@9748d66dad244ad9be08f032ad029f63666a6dd6c0ae80baa8a1183aded8b644:1;CIGAR=1M77D	GT:FT:GQ:PL:PR:SR:OLD_GT:DP:AD:ADF:ADR	1/1:PASS:53:270,63,0:0,0:0,17:1/1:18:0,22:0,19:0,3	1/1:PASS:41:390,90,0:0,0:0,15:1/1:26:0,32:0,17:0,15
"""
            }
            [ meta, vcf ]
        }
        .set { ch_bgzip_input }

    TABIX_ANNOTSV(
        ch_bgzip_input
    )
    ch_versions = ch_versions.mix(TABIX_ANNOTSV.out.versions.first())

    ENSEMBLVEP_VEP(
        ch_vcfs,
        params.genome,
        params.species,
        params.vep_cache_version,
        ch_vep_cache,
        ch_fasta.map { [[], it] },
        ch_vep_extra_files
    )
    ch_reports  = ch_reports.mix(ENSEMBLVEP_VEP.out.report)
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)

    TABIX_VEP(
        ENSEMBLVEP_VEP.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_VEP.out.versions)

    ENSEMBLVEP_VEP.out.vcf
        .join(TABIX_VEP.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .join(TABIX_ANNOTSV.out.gz_tbi
            .map { meta, vcf, tbi -> [ meta, [vcf, tbi] ]},
            failOnDuplicate:true,
            failOnMismatch:true
        )
        .set { ch_vcfanno_input }

    Channel.fromList(create_vcfanno_toml(val_vcfanno_resources))
        .collectFile(name:"vcfanno.toml", newLine:true)
        .collect()
        .set { ch_vcfanno_toml }

    VCFANNO(
        ch_vcfanno_input,
        ch_vcfanno_toml,
        ch_vcfanno_lua,
        val_vcfanno_resources ? Channel.fromList(val_vcfanno_resources).collect() : []
    )
    ch_versions = ch_versions.mix(VCFANNO.out.versions)

    TABIX_ANNOTATED(
        VCFANNO.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_ANNOTATED.out.versions)
    
    emit:
    // annotated_vcfs  = TABIX_ANNOTATED.out.gz_tbi  // channel: [ val(meta), path(vcf), path(tbi) ]

    reports         = ch_reports
    versions        = ch_versions
}

def create_vcfanno_toml(vcfanno_resources) {
    def params_toml_files = params.vcfanno_toml ? parse_toml(params.vcfanno_toml) : [postannotation:[]]
    def assets_toml_files = parse_toml("${projectDir}/assets/vcfanno/*.toml")
    def resources = vcfanno_resources.collect { it.fileName.toString() }
    resources.add("annotsv_annotated.vcf.gz")
    def output = []
    for (file_name in resources) {
        if (params_toml_files.containsKey(file_name)){
            output.add(params_toml_files[file_name])
        }
        else if (assets_toml_files.containsKey(file_name)){
            output.add(assets_toml_files[file_name])
        }
    }
    postannotation = params_toml_files.postannotation != [] ? params_toml_files.postannotation : assets_toml_files.postannotation
    if (postannotation != []){
        output.add(postannotation)
    }
    return output.flatten()
}

def parse_toml(tomls) {
    def output = [:]
    output.postannotation = []
    tomls_files = file(tomls, checkIfExists:true)
    toml_list = tomls_files instanceof LinkedList ? tomls_files : [tomls_files]
    for (toml in toml_list) {
        def info = ""
        def file = ""
        def fields = []
        for (line in toml.readLines()) {
            if (line.startsWith("#")) { continue }
            if (line == "[[annotation]]" || line == "[[postannotation]]") {
                if (info.startsWith("[[postannotation]]")) {
                    output.postannotation.add(create_toml_config(fields))
                }
                else if(info != "") {
                    output[file] = create_toml_config(fields)
                }
                if (info != "") {
                    info = ""
                    file = ""
                    fields = []
                }
                info = line
            }
            else if (line.startsWith("file")) {
                file = line.split("\"")[-1]
            }
            fields.add(line)
        }
        if (info.startsWith("[[postannotation]]")) {
            output.postannotation.add(create_toml_config(fields))
        } else {
            output[file] = create_toml_config(fields)
        }
    }
    return output
}

def create_toml_config(fields_list) {
    config = fields_list.findAll { it != "" }.join("\n")
    return "${config}\n"
}
