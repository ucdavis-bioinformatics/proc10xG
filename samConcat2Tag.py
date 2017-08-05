#!/usr/bin/env python
'''
Copyright 2017 Matt Settles
Created June 8, 2017

bwa mem -C option concatenats the fasta/fastq
CTACATTGTCAAGGGT:E00558:34:HGCJ3ALXX:1:1101:2108:1731   99      000000F 922571  60      127M    =       922961  517     ACTCGGGGAGGTGTTAGCTGCTGCCTCACACATTGGGTTTATAGGCTGAATCTTGTTCTCTTTAGGCTTCCAGAGTTTTCTCAGTTACTATTTCTCCTGTCACATACTCGCTGCTTCTTCTGTCATA JJJJJJ<JJF<7A7FJJJJJJ<JJJAJAJJJFJFFFJ----AJJFJ---7---<FJJ<JF<7FFFJJJFJJAJF-AAFFFFF-AFJF7FF<A--FJJJAF)-7-77<<7--)7)<<--77A7-<--< NM:i:3  MD:Z:74T34A3T13 AS:i:112        XS:i:19 1:N:0:GOOD:CCGATTAA:CTACATTGTCAAGGGT:<AAFFJJFJJFJJJJJ:CCAGTGA:J<FFFJJ

This pulls it out, 9 columns and produces new 10x tags in the bam then writes to out
'''
import sys
import os
from optparse import OptionParser  # http://docs.python.org/library/optparse.html

# TODO: ADD parameter for sample ID
usage = "usage: %prog -o output_base inputfile.SAM"
parser = OptionParser(usage=usage)
parser.add_option('-o', '--output_base', help="output file basename",
                  action="store", type="str", dest="output_base", default=None)

(options, args) = parser.parse_args()  # uncomment this line for command line support


if len(args) == 1:
    infile = args[0]
    # Start opening input/output files:
    if not os.path.exists(infile):
        sys.exit("Error, can't find input file %s" % infile)
    insam = open(infile, 'r')
else:
    # reading from stdin
    insam = sys.stdin

base = options.output_base

if base is None:
    out = sys.stdout
else:
    out = open(base + ".sam", 'w')

for line in insam:
    # Comment/header lines start with @
    if line[0] != "@" and len(line.strip().split()) > 2:
        line2 = line.strip().split()
        # Handle TAG:
        # get the final concatenatd tag
        tag = line2[-1]
        if (tag[0:6] in ['1:N:0:', '2:N:0:']):
            tsplit = tag.split(":")
            if len(tsplit) != 9:
                sys.exit("sam file has concatenated info, but its the wrong size")
            # fixed barcode
            line2 = line2[0:-1]
            line2.extend(["ST:Z:" + tsplit[3], "BX:Z:" + line2[0].split(":")[0] + "-1", "BC:Z:" + tsplit[4], "QT:Z:" + '!' * len(tsplit[4]), "RX:Z:" + tsplit[5], "QX:Z:" + tsplit[6], "TR:Z:" + tsplit[7], "TQ:Z:" + tsplit[8]])
            out.write('\t'.join(line2) + '\n')
        else:   # Does not contain a concatenated tag as expected by bwa mem
            out.write('\t'.join(line2) + '\n')
    else:  # Its the header lines, so just put back on the stream/file
        out.write(line)


if base is not None:
    out.close()
