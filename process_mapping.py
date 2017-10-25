#!/usr/bin/env python
"""
Copyright 2017 Matt Settles
Created June 8, 2017
"""
from optparse import OptionParser
import os
import sys
import time
import traceback
import signal

from subprocess import Popen
from subprocess import PIPE


#def seqToHash(seq):
#    encoding = {'a': 0, 'c': 1, 'g': 2, 't': 3, 'A': 0, 'C': 1, 'G': 2, 'T': 3}
#    result = 0
#    i = 0
#    while i < len(seq):
#        result += encoding.get(seq[i]) * 4**i
#        i += 1
#    return result


def sp_bwa_index(ref, overwrite=False):
    if os.path.isfile(ref):
        if os.path.isfile(ref + '.sa') and not overwrite:
            sys.stderr.write('Found existing bwo index for %s\n' % ref)
            return 0
        else:
            FNULL = open(os.devnull, 'w')
            call = 'bwa index'
            call = call + ' ' + ref
            sys.stdout.write(call + '\n')
            p = Popen(['bwa index', ref],
                      stdout=FNULL,
                      stderr=FNULL,
                      bufsize=-1,
                      preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
            p.communicate()
            if p.returncode:
                sys.stderr.write('Something in bwa index went wrong\n')
                raise
            # system call, check for return
            sys.stderr.write('Successfully indexed %s\n' % ref)
            return 0
    else:
        sys.stderr.write("%s Reference file not found\n" % ref)
        return 1
    sys.stderr.write('Something in bwa index went wrong\n')
    raise


def sp_bwa_map_bc(reads, ref, overwrite=False, sensitivity=0, procs=1):
    # build the call,
    if sp_bwa_index(ref, overwrite) != 0:
        sys.exit(1)

    call = 'bwa mem -p -a -t ' + procs
    p = Popen([call, ref],
              stdin=PIPE,
              stdout=PIPE,
              stderr=None,
              bufsize=-1,
              shell=True,
              executable='/bin/bash',
              preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    p.communicate(input=reads)
    if p.returncode:
        raise
    return p.stdout


def reverseComplement(s):
    """
    given a seqeucne with 'A', 'C', 'T', and 'G' return the reverse complement
    """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    letters = list(s)
    try:
        letters = [basecomplement[base] for base in letters]
    except:
        raise
    return ''.join(letters[::-1])


def reverse(s):
    """
    given a sequence return the reverse
    """
    letters = list(s)
    return ''.join(letters[::-1])


class bcProcessing:
    orig_bc_lines = []
    remap_reads = []
    remap_lines = []
    ok_bc_lines = []

    def __init__(self, minMQ=40, verbose=False):
        self.verbose = verbose
        self.minMQ = minMQ

    def addLine(self, line):
        self.orig_bc_lines.append(line)

    def clearbc(self):
        self.orig_bc_lines = []
        self.remap_lines = []
        self.ok_bc_lines = []

    def addRead(self, read):
        """
        Add a pair of reads to the output queue
        """
        self.R1.append(read[0])
        self.R2.append(read[3])
        self.mcount += 1

    def mapReads(self):
        """
        Write the paired reads in the queue to the output files
        """
        if (len(self.R1) == 0):
            pass
        else:
            if not self.isOpen:
                try:
                    if self.open() == 1:
                        sys.stderr.write('ERROR:[IlluminaFourReadOutput] ERROR Opening files for writing\n')
                        raise
                except Exception:
                    raise
            try:
                self.R1f.write('\n'.join(self.R1) + '\n')
                self.R4f.write('\n'.join(self.R2) + '\n')
            except Exception:
                sys.stderr.write('ERROR:[IlluminaFourReadOutput] Cannot write reads to file with prefix: %s\n' % self.output_prefix)
                raise
            self.R1 = []
            self.R2 = []
            self.B1 = []
            self.B2 = []
            self.close()

    def process(self, refDB):
        """
        process reads against a reference fasta file
        """
        try:
            bc = self.orig_bc_lines[0].split(":")[0]
            mapped_pairs_count = 0
            remapped_pairs_count = 0
            mapped_singles_count = 0
            remapped_singles_count = 0
            secondary_alignment = 0
            count = 0
            PE1 = {}
            PE2 = {}

            for line in self.orig_bc_lines:
                # 0x1 template having multiple segments in sequencing
                # 0x2 each segment properly aligned according to the aligner
                # 0x4 segment unmapped
                # 0x8 next segment in the template unmapped
                # 0x10 SEQ being reverse complemented
                # 0x20 SEQ of the next segment in the template being reversed
                # 0x40 the first segment in the template
                # 0x80 the last segment in the template
                # 0x100 secondary alignment
                # 0x200 not passing quality controls
                # 0x400 PCR or optical duplicate
                lbc = line.split(":")[0]
                if lbc != bc:
                    sys.err("Something went wrong, more than one barcode in process barcodes")

                count += 1
                line2 = line.strip().split()
                flag = int(line2[1])
                # Secondary alignment
                if (flag & 0x100):  # not sure what to do with secondary alignment yet, for now ignore
                    secondary_alignment += 1
                    continue

                mq = int(line2[4])
#                if mq < self.minMQ:
#                    # add to those that need to be remapped
#                    self.remap_lines.append(line)
                # Handle SE:
                # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1
                if not (flag & 0x1):  # SE READ, shouldn't see singles, maybe handle them later
                    # if not (flag & 0x4 and mq >= self.minMQ):  # MAPPED
                    #     self.ok_bc_lines.append(line)
                    #     # TODO: NEED to determine read cloud for read
                    #     mapped_singles_count += 1
                    # else:  # UNMAPPED or Poor minMQ, remap the read
                    #     if (flag & 0x10):  # reverse complement
                    #         line2[9] = reverseComplement(line2[9])
                    #         line2[10] = reverse(line2[10])
                    #     self.addRead(['\n'.join(['@' + line2[0] + ' 1:N:O', line2[9], '+', line2[10]])])
                    #     remapped_singles_count += 1
                    continue
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

        except (KeyboardInterrupt, SystemExit):
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
            return 1
        except:
            sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1


def main(insam, outsam, output_all, verbose):
    global file_path
    refDict = {}
    proc_bc = bcProcessing(verbose)

    #gbcDict = {}
    # Load the gem barcode dictionary with the whitelist
    #with open(os.path.join(file_path, 'barcodes/4M-with-alts-february-2016.txt'), 'r') as f:
    #    for bc_sequence in f:
    #        gbcDict[seqToHash(bc_sequence.strip())] = bc_sequence.strip()
    #    if verbose:
    #        sys.stderr("Finished reading in barcode whitelist")

    line_count = 0
    current_bc = None
    current_bc_count = 0
    bc_count = 0
    for line in insam:
        # Comment/header lines start with @
        if line[0] != "@" and len(line.strip().split()) > 2:
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
            # pass header directly to output
            outsam.write(line)
            if line[0:3] == "@SQ":
                # reference sequence id
                sp = line.split()
                refDict[sp[1][3:]] = int(sp[2][3:])

        if line_count % 100000 == 0 and line_count > 0 and verbose:
            print "Records processed: %s" % (line_count)


#####################################
# Parse options and setup #
usage = "usage %prog -o [output file prefix (path + name)] -(a) --quiet samfile"
usage += "%prog will process alignment file produced by processing_10xReads and do some stuff, assumes sam file is sorted by read name"
parser = OptionParser(usage=usage, version="%prog 0.0.1")

parser.add_option('-o', '--output', help="Directory + filename to output sam file, or stdout",
                  action="store", type="str", dest="outfile", default="stdout")

parser.add_option('-a', '--all', help="output all alignment, not just those with valid gem barcode (STATUS is UNKNOWN, or AMBIGUOUS)",
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
    outsam = sys.stdout
else:
    outsam = open(outfile, 'r')

output_all = options.output_all
verbose = options.verbose

# need to check, can write to output folder

# global variables
file_path = os.path.dirname(os.path.realpath(__file__))

stime = time.time()

main(insam, outsam, output_all, verbose)

sys.exit(0)
