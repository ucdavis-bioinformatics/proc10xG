#!/usr/bin/env python
"""
Copyright 2017 Matt Settles
Created June 8, 2017
"""
from optparse import OptionParser
from collectons import Counter
import os
import sys
import time
import traceback
import signal

from subprocess import Popen
from subprocess import PIPE

                # Handle PE:
                # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
                if (flag & 0x1):  # PE READ
                    if (not (flag & 0x4) and not (flag & 0x8)):  # both pairs mapped
                        if (flag & 0x40):  # is this PE1 (first segment in template)
                            # PE1 read, check that PE2 is in dict
                            ID = line2[0]
                            if ID in PE2:
                                if mq >= self.minMQ and int(PE2[ID].strip().split()[4]) >= self.minMQ:  # check MQ of both reads
                                    self.ok_bc_lines.append(line)
                                    self.ok_bc_lines.append(PE2[ID])
                                    del PE2[ID]
                                    # TODO: NEED to determine read cloud for read
                                    mapped_pairs_count += 1
                                else:
                                    if (flag & 0x10):  # reverse complement
                                        line2[9] = reverseComplement(line2[9])
                                        line2[10] = reverse(line2[10])
                                    r1 = '\n'.join(['@' + line2[0] + ' 1:N:O', line2[9], '+', line2[10]])  # sequence + qual
                                    rl2 = PE2[ID].strip().split()
                                    if (int(rl2[1]) & 0x10):  # reverse complement
                                        rl2[9] = reverseComplement(rl2[9])
                                        rl2[10] = reverse(rl2[10])
                                    r2 = '\n'.join(['@' + rl2[0] + ' 2:N:O', rl2[9], '+', rl2[10]])  # sequence + qual
                                    self.addRead('\n'.join([r1, r2]))
                                    del PE2[ID]
                                    remapped_pairs_count += 1
                            else:
                                PE1[ID] = line
                        elif (flag & 0x80):  # is this PE2 (last segment in template)
                            # PE2 read, check that PE1 is in dict and write out
                            ID = line2[0]
                            if ID in PE1:
                                if mq >= self.minMQ and int(PE1[ID].strip().split()[4]) >= self.minMQ:  # check MQ of both reads
                                    self.ok_bc_lines.append(line)
                                    self.ok_bc_lines.append(PE1[ID])
                                    del PE1[ID]
                                    # TODO: NEED to determine read cloud for read
                                    mapped_pairs_count += 1
                                else:
                                    if (flag & 0x10):  # reverse complement
                                        line2[9] = reverseComplement(line2[9])
                                        line2[10] = reverse(line2[10])
                                    r2 = '\n'.join(['@' + line2[0] + ' 2:N:O', line2[9], '+', line2[10]])  # sequence + qual
                                    rl1 = PE1[ID].strip().split()
                                    if (int(rl1[1]) & 0x10):  # reverse complement
                                        rl1[9] = reverseComplement(rl1[9])
                                        rl1[10] = reverse(rl1[10])
                                    r1 = '\n'.join(['@' + rl1[0] + ' 1:N:O', rl1[9], '+', rl1[10]])  # sequence + qual
                                    self.addRead('\n'.join([r1, r2]))
                                    del PE1[ID]
                                    remapped_pairs_count += 1
                            else:
                                PE2[ID] = line
                    else:  # an 'unmapped' pair, at least 1 unmapped
                        if (flag & 0x40):  # is this PE1 (first segment in template)
                            # PE1 read, check that PE2 is in dict and write out
                            ID = line2[0]
                            if ID in PE2:
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r1 = '\n'.join(['@' + line2[0] + ' 1:N:O', line2[9], '+', line2[10]])  # sequence + qual
                                rl2 = PE2[ID].strip().split()
                                if (int(rl2[1]) & 0x10):  # reverse complement
                                    rl2[9] = reverseComplement(rl2[9])
                                    rl2[10] = reverse(rl2[10])
                                r2 = '\n'.join(['@' + rl2[0] + ' 2:N:O', rl2[9], '+', rl2[10]])  # sequence + qual
                                self.addRead('\n'.join([r1, r2]))
                                del PE2[ID]
                                remapped_pairs_count += 1
                            else:
                                PE1[ID] = line
                        elif (flag & 0x80):  # is this PE2 (last segment in template)
                            # PE2 read, check that PE1 is in dict and write out
                            ID = line2[0]
                            if ID in PE1:
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r1 = '\n'.join(['@' + line2[0] + ' 1:N:O', line2[9], '+', line2[10]])  # sequence + qual
                                rl2 = PE2[ID].strip().split()
                                if (int(rl2[1]) & 0x10):  # reverse complement
                                    rl2[9] = reverseComplement(rl2[9])
                                    rl2[10] = reverse(rl2[10])
                                r2 = '\n'.join(['@' + rl2[0] + ' 2:N:O', rl2[9], '+', rl2[10]])  # sequence + qual
                                self.addRead('\n'.join([r1, r2]))
                                del PE2[ID]
                                remapped_pairs_count += 1
                            else:
                                PE2[ID] = line


def main(insam, output_all, verbose):
    global file_path
    refDict = {}
    bcDict = {}

    line_count = 0
    bc_count = 0

    for line in insam:
        # Comment/header lines start with @
        if line[0] == "@":
            # pass header directly to output
            if line[0:3] == "@SQ":
                # reference sequence id
                sp = line.split()
                refDict[sp[1][3:]] = int(sp[2][3:])
        elif line[0] != "@" and len(line.strip().split()) > 2:
            line_count += 1
            bc = line.split(":")[0]

            # instead check the ST:Z:GOOD for GOOD or MATCH or MISMATCH1
            if line.split()[15][5:] not in ['GOOD', 'MATCH', 'MISMATCH1']:
                # if seqToHash(bc) not in gbcDict:
                # barcode does not match whitelist
                if output_all:
                    # if output_all pass line directly to output
                    outsam.write(line)
            elif bc == current_bc:
                # add line to bc processing
                proc_bc.addLine(line)
                current_bc_count += 1
            elif current_bc is None:
                current_bc = bc
                # add line to bc processing
                proc_bc.addLine(line)
                current_bc_count += 1
            else:
                # this is a new barcode
                # can add a check to see if seen bc before, which is a no-no
                # process the bc
                proc_bc.process()
                # record having processed the barcode
                # output to sam file
                bc_count += 1
                proc_bc.clearbc()
                current_bc = bc
                # add line to bc processing
                current_bc_count = 1
                proc_bc.addLine(line)
        else:
            # Not sure what happened
            sys.stderr.write("Unknown line: %s" % line)

        if line_count % 100000 == 0 and line_count > 0 and verbose:
            print "Records processed: %s" % (line_count)


#####################################
# Parse options and setup #
usage = "usage %prog -o [output file prefix (path + name)] -(a) --quiet samfile"
usage += "%prog will process alignment file produced by processing_10xReads and do profile each barcode"
parser = OptionParser(usage=usage, version="%prog 0.0.1")

parser.add_option('-o', '--output', help="Directory + filename to output bc stats",
                  action="store", type="str", dest="outfile", default="bc_profile.txt")

parser.add_option('-a', '--all', help="output all barcodes, not just those with valid gem barcode (STATUS is UNKNOWN, or AMBIGUOUS)",
                  action="store_true", dest="output_all", default=False)

parser.add_option('--quiet', help="turn off verbose output",
                  action="store_false", dest="verbose", default=True)

(options, args) = parser.parse_args()

if len(args) == 1:
    infile = args[0]
    # Start opening input/output files:
    if not os.path.exists(infile):
        sys.exit("Error, can't find input file %s" % infile)
    insam = open(infile, 'r')
else:
    # reading from stdin
    insam = sys.stdin


outfile = options.outfile
if outfile == "stdout":
    outf = sys.stdout
else:
    outf = open(outfile, 'r')

output_all = options.output_all
verbose = options.verbose

# need to check, can write to output folder

# global variables
file_path = os.path.dirname(os.path.realpath(__file__))

stime = time.time()

main(insam, outf, output_all, verbose)

sys.exit(0)
