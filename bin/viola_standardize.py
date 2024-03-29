#!/usr/local/bin/python

import argparse
import os
import shutil
import numpy as np
import pandas as pd

import viola

if __name__ == "__main__":
    # Setting up argparser
    parser = argparse.ArgumentParser(description="A script to standardize VCFs using Viola-SV")
    parser.add_argument("vcf", metavar="FILE", type=str, help="The called VCF")
    parser.add_argument("caller", metavar="STRING", type=str, help="The caller used to call the VCF")
    parser.add_argument("out_file", metavar="FILE", type=str, help="The standardized VCF")
    parser.add_argument("patient_name", metavar="STRING", type=str, help="The name of the patient in the VCF file")

    args = parser.parse_args()

    in_file = args.vcf
    caller = args.caller
    out_file = args.out_file
    patient_name = args.patient_name

    if caller == "smoove":
        caller = "lumpy"

    if caller == "gridss":
        svlen_not_added = True
        old_vcf = f"old_{in_file}"
        os.rename(in_file, old_vcf)
        with open(old_vcf, "r") as old:
            with open(in_file, "w") as new:
                for line in old.readlines():
                    if line.startswith("##INFO") and svlen_not_added:
                        svlen_not_added = False
                        new.write(
                            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="The length of the structural variant.">\n'
                        )
                    new.write(line.replace("CIRPOS", "CIEND"))
        raw_vcf = viola.read_vcf(in_file, variant_caller=caller, patient_name=patient_name)
        vcf = raw_vcf.breakend2breakpoint()

        # Genotype the variants
        # TODO Fix the calculation of the genotypes to be more precise
        # TODO Analyse the amount of variants are inside the intervals 0-20 >20-40 >40-60 >60-80 >80
        genotype = vcf.get_table("formats")[["id", "sample", "format", "value"]].rename(columns={"value": "af"})
        genotypes = genotype[genotype["format"] == "AF"]
        genotypes.loc[:, "gt"] = np.where(
            genotypes.loc[:, "af"] >= 0.75, "1/1", np.where(genotypes.loc[:, "af"] <= 0.25, "0/0", "0/1")
        )
        genotypes.loc[:, "format"] = "GT"
        formats = pd.merge(vcf.get_table("formats"), genotypes, how="left", on=["format", "id", "sample"])
        formats.loc[formats["gt"].notna(), "value"] = formats["gt"]
        formats_done = formats.drop(columns=["gt", "af"])
        vcf.replace_table("formats", formats_done)

        # Fix the ALT fields for breakpoint notation
        orig_table = vcf.get_table("positions")
        alt_table = viola.io.parser.create_alt_field_from_position(orig_table)
        alt_table["alt"] = np.where(alt_table["svtype"] == "BND", orig_table["alt"], alt_table["alt"])
        vcf.replace_table("positions", alt_table)
    else:
        vcf = viola.read_vcf(in_file, variant_caller=caller, patient_name=patient_name)

    # Fix IDs to contain the caller
    ids = vcf.ids
    new_ids = []
    for id in ids:
        number = id.split(":")[-1]
        new_ids.append(f"{caller}_{number}")
    vcf.replace_svid(ids, new_ids)

    # Write to output file
    try:
        vcf.to_vcf(out_file)
    except TypeError as e:
        print("No variants found, returning a copy of the input file")
        shutil.copyfile(in_file, out_file)
