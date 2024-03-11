#!/usr/bin/env python

import gzip
import argparse

def main() -> None:
    parser = argparse.ArgumentParser(description="Extract the callers from a VCF file and update the CALLERS header field")
    parser.add_argument('input', metavar='FILE', type=str, help="The VCF file to process")
    parser.add_argument('output', metavar='FILE', type=str, help="The output VCF file with updated CALLERS header field")

    args = parser.parse_args()
    input = args.input

    if input.endswith(".gz"):
        with gzip.open(input, "rt") as file:
            process_file(file, args.output)
    else:
        with open(input, 'r') as file:
            process_file(file, args.output)
                
def process_file(file, output) -> None:
    with open(output, 'w') as out_file:
        for line in file:
            if line.startswith("#"):
                out_file.write(line)
            else:
                pre_info = "\t".join(line.split("\t")[0:7])
                info = line.split("\t")[7]
                after_info = "\t".join(line.split("\t")[8:])
                info_dict = {}
                for field in info.split(";"):
                    split_field = field.split("=", 1)
                    key, value = (split_field[0], split_field[1]) if len(split_field) == 2 else (split_field[0], "")
                    print(key, value)
                    info_dict[key] = value
                
                info_dict["CALLERS"] = ",".join(set([it.split("_")[0] for it in info_dict["IDLIST_EXT"].replace(".", ",").split(",")]))
                if len(info_dict["CALLERS"]) == 0:
                    del info_dict["CALLERS"]
                out_file.write(f"{pre_info}\t{';'.join([f'{key}={value}' for key, value in info_dict.items()])}\t{after_info}")

if __name__ == "__main__":
    main()

