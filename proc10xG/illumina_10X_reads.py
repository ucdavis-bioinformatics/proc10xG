# illumina_10X_reads.py
import sys
import os
import glob

from .misc import sp_gzip_read

class TwoRead10xLLRun:
    """
    Class to open/close and read a two read illumina 10X Genomics linked reads
    sequencing run. Data is expected to be in fastq format (possibly gzipped)
    """
    def __init__(self, read1, read2, gbctrim=16, prtrim=7, verbose=True):
        """
        Initialize a TwoRead10xRun object with expandible paths (with glob) to the two
        sequencing read files. A vector of multiple files per read is allowed.
        """
        self.verbose = verbose
        self.gbctrim = gbctrim
        self.trim = prtrim
        self.mcount = 0
        self.fread1 = []
        self.fread2 = []
        try:
            for fread in read1:
                self.fread1.extend(glob.glob(fread))
                if len(self.fread1) == 0 or not all(os.path.isfile(f) for f in self.fread1):
                    sys.stderr.write('ERROR[TwoRead10xLLRun]\tread1 file(s) not found\n')
                    raise Exception
            for fread in read2:
                self.fread2.extend(glob.glob(fread))
                if len(self.fread2) == 0 or not all(os.path.isfile(f) for f in self.fread2):
                    sys.stderr.write('ERROR[TwoRead10xLLRun]\tread2 file not found\n')
                    raise Exception
            if len(self.fread1) != len(self.fread2):
                sys.stderr.write('ERROR[TwoRead10xLLRun]\tInconsistent number of files for each read\n')
                raise Exception
        except Exception:
            raise
        # record the number of files per read
        self.numberoffiles = len(self.fread1)
        self.isOpen = False

    def open(self):
        """
        Open a OneReadIlluminaRun file set, if file ends in .gz, open will use gzip
        """
        if self.isOpen:
            self.close()
        if self.numberoffiles > 0:
            try:
                read1 = self.fread1.pop()
                if read1.split(".")[-1] == "gz":
                    self.R1 = sp_gzip_read(read1)
                else:
                    self.R1 = open(read1, 'r')
                read2 = self.fread2.pop()
                if read2.split(".")[-1] == "gz":
                    self.R2 = sp_gzip_read(read2)
                else:
                    self.R2 = open(read2, 'r')
            except Exception:
                sys.stderr.write('ERROR[TwoRead10xLLRun]\tcannot open input files\n')
                raise
            self.isOpen = True
            self.numberoffiles -= 1
            if self.verbose:
                sys.stderr.write("PROCESS\tFILES\t{0},{1}\n".format(read1, read2))
            return 0
        else:
            return 1

    def close(self):
        """
        Close a TwoRead10xRun file set
        """
        self.R1.close()
        self.R2.close()
        self.isOpen = False

    def count(self):
        """
        Provide the current count of reads read
        """
        return self.mcount

    def nfiles(self):
        """
        provide the number of files left to process (unopened)
        """
        return self.numberoffiles

    def next_raw(self, whitelist, ncount=1):
        """
        Extract and store the next [count] reads into a TwoSequenceReadSet object.
        If the file object is not open, or if 'next' reaches the end of a file, it will
        attempt to open the file in the list, or gracefully exit
        """
        if not self.isOpen:
            try:
                if self.open() == 1:
                    sys.stderr.write('PROCESS\tERROR:[TwoRead10xLLRun] ERROR Opening files for reading\n')
                    raise
            except Exception:
                raise
        reads = []
        i = 0
        while i < ncount:
            try:
                # pull in read 1
                status = 'UNKNOWN'
                id1 = next(self.R1).strip()
                seq1 = next(self.R1).strip()
                next(self.R1)  # *
                qual1 = next(self.R1).strip()
                assert(len(seq1) == len(qual1))
                if id1 == '' or seq1 == ''or qual1 == '':
                    self.close()
                    raise StopIteration
                # pull in read2
                id2 = next(self.R2).strip()
                seq2 = next(self.R2).strip()
                next(self.R2)  # *
                qual2 = next(self.R2).strip()
                assert(len(seq2) == len(qual2))
                if id2 == '' or seq2 == ''or qual2 == '':
                    self.close()
                    raise StopIteration
                # check to make sure the IDs match across the two reads
                assert(id1.split()[0] == id2.split()[0])
                rid = id1.split()[0][1:]
                rbc = (id1.split()[1]).split(':')[3]
                if rbc == '':
                    rbc = "0"
                gbc = seq1[0:self.gbctrim]
                gbcq = qual1[0:self.gbctrim]
                trim = seq1[self.gbctrim:self.gbctrim + self.trim]
                trimq = qual1[self.gbctrim:self.gbctrim + self.trim]
                seq1 = seq1[self.gbctrim + self.trim:]
                qual1 = qual1[self.gbctrim + self.trim:]
                # get the whitelisted gem barcode sequence, (None) if UNKNOWN
                # or AMBIGUOUS, and status (MATCH, MISMATCH1, AMBIGUOUS, UNKNOWN)
                gem_bc, status = whitelist.check_barcode(gbc)

                fragment = {'id': rid,
                            'status': status,
                            'gem_bc': gem_bc,
                            'sgem_bc': gbc,
                            'sgem_qual': gbcq,
                            'trim_seq': trim,
                            'trim_qual': trimq,
                            'read1_seq': seq1,
                            'read1_qual': qual1,
                            'read2_seq': seq2,
                            'read2_qual': qual2,
                            'library_bc': rbc}
                reads.append(fragment)
                self.mcount += 1
            except StopIteration:
                if self.numberoffiles > 0:
                    try:
                        if self.open() == 1:
                            sys.stderr.write('ERROR[TwoRead10xLLRun]\terror opening files for reading\n')
                            raise
                    except Exception:
                        raise Exception
                    continue
                raise StopIteration
            except Exception:
                sys.stderr.write('ERROR[TwoRead10xLLRun]\terror reading next read\n')
                raise
            i += 1
        return reads


class TwoRead10xLLOutput:
    """
    Given output_prefix filename and reads, output them to paired files (possibly gzipped)
    """
    def __init__(self, output_prefix, uncompressed=False, interleaved=False, verbose=True):
        """
        Initialize an TwoRead10xOutput object with output_prefix and whether or not
        output should be compressed with gzip [uncompressed True/False]
        """
        self.isOpen = False
        self.output_prefix = output_prefix
        self.interleaved = interleaved
        self.uncompressed = uncompressed
        self.mcount = 0
        self.verbose = True

        if output_prefix == "stdout":
            self.interleaved = True
            self.uncompressed = True
        elif self.uncompressed is True:
            if os.path.isfile(self.output_prefix + "_R1_001.fastq"):
                if self.verbose:
                    sys.stderr.write('PROCESS\tWARNING[TwoRead10xLLOutput]\tFile with prefix: {0} exists, DELETING\n'.format(self.output_prefix + "_R1_001.fastq"))
                try:
                    if self.interleaved:
                        os.remove(self.output_prefix + "_R1_001.fastq")
                    else:
                        os.remove(self.output_prefix + "_R1_001.fastq")
                        os.remove(self.output_prefix + "_R2_001.fastq")
                except Exception:
                    sys.stderr.write('ERROR[TwoRead10xLLOutput]\tCannot delete file with prefix: {0}\n'.format(self.output_prefix + "_R1_001.fastq"))
                    raise
        else:
            if os.path.isfile(self.output_prefix + "_R1_001.fastq.gz"):
                if self.verbose:
                    sys.stderr.write('PROCESS\tWARNING[TwoRead10xLLOutput]\tFile with prefix: {0} exists, DELETING\n'.format(self.output_prefix + "_R1_001.fastq"))
                try:
                    if self.interleaved:
                        os.remove(self.output_prefix + "_R1_001.fastq.gz")
                    else:
                        os.remove(self.output_prefix + "_R1_001.fastq.gz")
                        os.remove(self.output_prefix + "_R2_001.fastq.gz")
                except Exception:
                    sys.stderr.write('ERROR[TwoRead10xLLOutput] Cannot delete file with prefix: {0}\n'.format(self.output_prefix + "_R1_001.fastq"))
                    raise

    def open(self):
        """
        Open the two read files for writing, appending _R1.fastq and _R2.fastq to the output_prefix.
        Create directories as needed.
        """
        if self.isOpen:
            self.close()
        try:
            if self.output_prefix == "stdout":
                self.interleaved = True
                self.R1f = sys.stdout
            else:
                make_sure_path_exists(os.path.dirname(self.output_prefix))
                if self.uncompressed is True:
                    self.R1f = open(self.output_prefix + '_R1_001.fastq', 'w')
                    if not self.interleaved:
                        self.R2f = open(self.output_prefix + '_R2_001.fastq', 'w')
                else:
                    self.R1f = sp_gzip_write(self.output_prefix + '_R1_001.fastq.gz')
                    if not self.interleaved:
                        self.R2f = sp_gzip_write(self.output_prefix + '_R2_001.fastq.gz')
        except Exception:
            sys.stderr.write('ERROR[TwoRead10xLLOutput] Cannot write reads to file with prefix: {0}\n'.format(self.output_prefix))
            raise
        self.isOpen = True
        return 0

    def close(self):
        """
        Close an TwoRead10xOutput file set
        """
        try:
            self.R1f.close()
            if not self.interleaved:
                self.R2f.close()
        except Exception:
            raise
        self.isOpen = False
        if self.verbose:
            sys.stderr.write("PROCESS\tFILES\tWrote {0} reads to output\n".format(self.mcount))

    def count(self):
        """
        Provide the current read count for the file output
        """
        return self.mcount

    def writeProcessedPairedFastq(self, fragment):
        newid = '@' + (':').join([fragment['gem_bc'], fragment['id'], fragment['sgem_bc'], fragment['sgem_qual'], fragment['trim_seq'], fragment['trim_qual'], fragment['status']])
        # read 1
        self.R1f.write((' ').join([newid, (':').join(['1', 'N', '0', fragment['library_bc']])]) + '\n')
        self.R1f.write(fragment['read1_seq'] + '\n')
        self.R1f.write('+\n')
        self.R1f.write(fragment['read1_qual'] + '\n')
        # read 2
        self.R2f.write((' ').join([newid, (':').join(['2', 'N', '0', fragment['library_bc']])]) + '\n')
        self.R2f.write(fragment['read2_seq'] + '\n')
        self.R2f.write('+\n')
        self.R2f.write(fragment['read2_qual'] + '\n')
        self.mcount += 1

    def writeProcessedFastqInterleaved(self, fragment):
        newid = '@' + (':').join([fragment['gem_bc'], fragment['id'], fragment['sgem_bc'], fragment['sgem_qual'], fragment['trim_seq'], fragment['trim_qual'], fragment['status']])
        # read 1
        self.R1f.write((' ').join([newid, (':').join(['1', 'N', '0', fragment['library_bc']])]) + '\n')
        self.R1f.write(fragment['read1_seq'] + '\n')
        self.R1f.write('+\n')
        self.R1f.write(fragment['read1_qual'] + '\n')
        # read 2
        self.R1f.write((' ').join([newid, (':').join(['2', 'N', '0', fragment['library_bc']])]) + '\n')
        self.R1f.write(fragment['read2_seq'] + '\n')
        self.R1f.write('+\n')
        self.R1f.write(fragment['read2_qual'] + '\n')

    def writeProcessedRead(self, fragments):
        """
        Write the paired read in the queue to the output files
        """
        if (len(fragments) == 0):
            pass
        else:
            if not self.isOpen:
                try:
                    if self.open() == 1:
                        sys.stderr.write('ERROR[TwoRead10xLLOutput]\tERROR Opening files for writing\n')
                        raise
                except Exception:
                    raise
            try:
                for frag in fragments:
                    if self.interleaved:
                        self.writeProcessedFastqInterleaved(frag)
                    else:
                        self.writeProcessedPairedFastq(frag)
            except Exception:
                sys.stderr.write('ERROR[TwoRead10xLLOutput]\tCannot write reads to file with prefix: {0}\n'.format(self.output_prefix))
                raise
