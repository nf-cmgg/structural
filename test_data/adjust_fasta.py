# A little script used to format the FASTA file correctly

from math import ceil
import re

asked_length = 61

with open("genome.fasta", "r") as input:
    with open("genome_adjusted.fasta", "w") as output:
        for line in input:
            length = len(line)
            if length > asked_length:
                lines_to_print = ceil(length / asked_length)
                for i in range(0,lines_to_print):
                    output.write(f'{line[i*(asked_length-1):(i+1)*(asked_length-1)].upper()}\n')
            elif re.match("^>chr\d\d?", line): 
                output.write(f'{line}')
            else:
                output.write(f'{line.upper()}')