#!/usr/local/bin/python

import argparse
import os

import viola

if __name__ == "__main__":
    # Setting up argparser
    parser = argparse.ArgumentParser(description="A script to standardize VCFs using Viola-SV")
    parser.add_argument('vcf', metavar='FILE', type=str, help="The called VCF")
    parser.add_argument('caller', metavar='STRING', type=str, help="The caller used to call the VCF")
    parser.add_argument('out_file', metavar='FILE', type=str, help="The standardized VCF")
    parser.add_argument('patient_name', metavar='STRING', type=str, help="The name of the patient in the VCF file")

    args = parser.parse_args()

    vcf = args.vcf
    caller = args.caller
    out_file = args.out_file
    patient_name = args.patient_name

    if caller == "smoove": caller = "lumpy"

    if caller == "gridss":
        svlen_not_added = True
        old_vcf = f'old_{vcf}'
        os.rename(vcf, old_vcf)
        with open(old_vcf, 'r') as old:
            with open(vcf, 'w') as new:
                for line in old.readlines():
                    if line.startswith("##INFO") and svlen_not_added:
                        svlen_not_added = False
                        new.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"The length of the structural variant.\">\n")
                    new.write(line.replace("CIRPOS", "CIEND"))

    viola.read_vcf(vcf, variant_caller=caller, patient_name=patient_name).breakend2breakpoint().to_vcf(out_file)
