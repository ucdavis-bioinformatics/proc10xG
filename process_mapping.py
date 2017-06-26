
import os
import sys
import time
import traceback
import signal

from subprocess import Popen
from subprocess import PIPE

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


class screeningApp:

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, fastq_file3, reference, overwrite, sensitivity, output_prefix, strict, procs, uncompressed=False, verbose=True, debug=False):
        """
            screen reads against a reference fasta file
        """
        self.verbose = verbose
        try:
            mapped_pairs_count = 0
            unmapped_pairs_count = 0
            mapped_singles_count = 0
            unmapped_singles_count = 0
            secondary_alignment = 0

            # Set up output
            self.run_out = {}
            self.run_out["mapped_pairs"] = IlluminaTwoReadOutput(output_prefix + '.mapped', uncompressed)
            self.run_out["unmapped_pairs"] = IlluminaTwoReadOutput(output_prefix + '.unmapped', uncompressed)
            self.run_out["mapped_singles"] = IlluminaOneReadOutput(output_prefix + '.mapped', uncompressed)
            self.run_out["unmapped_singles"] = IlluminaOneReadOutput(output_prefix + '.unmapped', uncompressed)

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
            PE1 = {}
            PE2 = {}

            i = 0
            for line in sp_bowtie2_screen(fastq_file1, fastq_file2, fastq_file3, reference, overwrite, sensitivity, procs):
                if i % 100000 == 0 and i > 0:
                    for key in self.run_out:
                        self.run_out[key].writeReads()
                    if self.verbose:
                        sys.stdout.write("Processed: %s, PE in ref: %s, SE in ref: %s\n" % (i, mapped_pairs_count, mapped_singles_count))
                if line[0] != "@":  # header line
                    i += 1

                    line2 = line.strip().split()
                    flag = int(line2[1])
                    # Secondary alignment
                    if (flag & 0x100):
                        secondary_alignment += 1
                        continue

                    # Handle SE:
                    # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1
                    if not (flag & 0x1):  # SE READ
                        if not (flag & 0x4):  # MAPPED
                            ID = line2[0]
                            if (flag & 0x10):  # reverse complement
                                line2[9] = reverseComplement(line2[9])
                                line2[10] = reverse(line2[10])
                            self.run_out["mapped_singles"].addRead(['\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])])
                            mapped_singles_count += 1
                        else:  # UNMAPPED
                            ID = line2[0]
                            if (flag & 0x10):  # reverse complement
                                line2[9] = reverseComplement(line2[9])
                                line2[10] = reverse(line2[10])
                            self.run_out["unmapped_singles"].addRead(['\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])])
                            unmapped_singles_count += 1
                        continue
                    # Handle PE:
                    # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
                    if (flag & 0x1):  # PE READ
                        if ((strict and not (flag & 0x4) and not (flag & 0x8)) or  # both pairs mapped
                                (not strict and (not (flag & 0x4) or not (flag & 0x8)))):  # at least one of the pair mapped
                            if (flag & 0x40):  # is this PE1 (first segment in template)
                                # PE1 read, check that PE2 is in dict and write out
                                ID = line2[0].split('_:_')[0]
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r1 = '\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])  # sequence + qual
                                if ID in PE2:
                                    self.run_out["mapped_pairs"].addRead([r1, PE2[ID]])
                                    del PE2[ID]
                                    mapped_pairs_count += 1
                                else:
                                    PE1[ID] = r1
                            elif (flag & 0x80):  # is this PE2 (last segment in template)
                                # PE2 read, check that PE1 is in dict and write out
                                ID = line2[0].split('_:_')[0]
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r2 = '\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])
                                if ID in PE1:
                                    self.run_out["mapped_pairs"].addRead([PE1[ID], r2])
                                    del PE1[ID]
                                    mapped_pairs_count += 1
                                else:
                                    PE2[ID] = r2
                        else:  # an 'unmapped' pair
                            if (flag & 0x40):  # is this PE1 (first segment in template)
                                # PE1 read, check that PE2 is in dict and write out
                                ID = line2[0].split('_:_')[0]
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r1 = '\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])  # sequence + qual
                                if ID in PE2:
                                    self.run_out["unmapped_pairs"].addRead([r1, PE2[ID]])
                                    del PE2[ID]
                                    unmapped_pairs_count += 1
                                else:
                                    PE1[ID] = r1
                            elif (flag & 0x80):  # is this PE2 (last segment in template)
                                # PE2 read, check that PE1 is in dict and write out
                                ID = line2[0].split('_:_')[0]
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r2 = '\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])
                                if ID in PE1:
                                    self.run_out["unmapped_pairs"].addRead([PE1[ID], r1])
                                    del PE1[ID]
                                    unmapped_pairs_count += 1
                                else:
                                    PE2[ID] = r2
            # Write out reads
            for key in self.run_out:
                self.run_out[key].writeReads()

            sys.stdout.write("Records processed: %s, PE in ref: %s, SE in ref: %s\n" % (i, mapped_pairs_count, mapped_singles_count))

            self.clean()
            return 0

        except (KeyboardInterrupt, SystemExit):
            self.clean()
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
            return 1
        except:
            self.clean()
            if not debug:
                sys.stderr.write("A fatal error was encountered. trying turning on debug\n")
            if debug:
                sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1

    def clean(self):
        if self.verbose:
            sys.stderr.write("Cleaning up.\n")
        try:
            pass
        except:
            pass


def main(samfile, output_dir, output_all, verbose):
    line_count = 0
    for line in samfile:
        if line_count % 100000 == 0 and line_count > 0 and verbose:
            print "Records processed: %s" % (line_count)
        # Comment/header lines start with @
        if line[0] != "@" and len(line.strip().split()) > 2:
            line_count += 1
            line2 = line.strip().split()
            flag = int(line2[1])
            # Handle SE:
            if (flag & 0x100):  # secondary alignment
                continue
            # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1
            if not (flag & 0x1) and not (flag & 0x4):
                ID = line2[0].split("#")[0]
                if (flag & 0x10):
                    line2[9] = reverseComplement(line2[9])
                    line2[10] = reverse(line2[10])
                outSE.write("@" + ID + '\n')
                outSE.write(line2[9] + '\n')
                outSE.write('+\n' + line2[10] + '\n')
                SE_written += 1
                continue
            # Handle PE:
            #logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
            if ((strict and (flag & 0x1) and not (flag & 0x4) and not (flag & 0x8))
                    or (not strict and (flag & 0x1) and (not (flag & 0x4) or not (flag & 0x8)))):
                if (flag & 0x40):  # is this PE1 (first segment in template)
                    #PE1 read, check that PE2 is in dict and write out
                    ID = line2[0].split("#")[0]
                    if (flag & 0x10):
                        line2[9] = reverseComplement(line2[9])
                        line2[10] = reverse(line2[10])
                    r1 = [line2[9], line2[10]]  # sequence + qual
                    if ID in PE2:
                        writeread(ID, r1, PE2[ID])
                        del PE2[ID]
                        PE_written += 1
                    else:
                        PE1[ID] = r1
                elif (flag & 0x80):  # is this PE2 (last segment in template)
                    #PE2 read, check that PE1 is in dict and write out
                    ID = line2[0].split("#")[0]
                    if (flag & 0x10):
                        line2[9] = reverseComplement(line2[9])
                        line2[10] = reverse(line2[10])
                    r2 = [line2[9], line2[10]]
                    if ID in PE1:
                        writeread(ID, PE1[ID], r2)
                        del PE1[ID]
                        PE_written += 1
                    else:
                        PE2[ID] = r2
        else:
            # pass header directly to output
            print "hello"


#####################################
# Parse options and setup #
usage = "usage %prog -o [output file prefix (path + name)] -(pbtg) --quiet samfile"
usage += "%prog will process alignment file produced by processing_10xReads and do some stuff"
parser = OptionParser(usage=usage, version="%prog 0.0.1")

parser.add_option('-o', '--output', help="Directory + filename to output sam file, or stdout",
                  action="store", type="str", dest="outfile", default="stdout")

parser.add_option('-a', '--all', help="output all alignment, not just those with valid gem barcode, STATUS will be UNKNOWN, or AMBIGUOUS",
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
    outfile = sys.stdout
else:
    outsam = open(outfile, 'r')

output_all = options.output_all
verbose = options.verbose

# need to check, can write to output folder

# global variables
file_path = os.path.dirname(os.path.abspath(__file__))

stime = time.time()

main(infile, outsam, output_all, verbose)

sys.exit(0)
