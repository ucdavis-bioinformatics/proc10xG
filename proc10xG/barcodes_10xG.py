# barcodes_10xG.py
# 10x Gem Barcode processing

import sys
from collections import Counter


class barcodes_10xG:
    __basedict = {
        'A': ['C', 'G', 'T'],
        'C': ['A', 'G', 'T'],
        'G': ['A', 'C', 'T'],
        'T': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T'],
        'a': ['C', 'G', 'T'],
        'c': ['A', 'G', 'T'],
        'g': ['A', 'C', 'T'],
        't': ['A', 'C', 'G'],
        'n': ['A', 'C', 'G', 'T']}
    """
    Class to store and then search and report on a 10X barcode whitelist
    """
    def __init__(self, whitelist_file, bc_trim, verbose):
        """
        Initialize the class with the barcode whitelist file
        """
        self.whitelist = whitelist_file
        self.bc_length = bc_trim
        self.verbose = verbose
        self._gbcDict = set()
        self._gbcCounter = Counter()
        self.statusCounter = {"MATCH": 0, "MISMATCH1": 0, "AMBIGUOUS": 0, "UNKNOWN": 0}

        # Load the gem barcode dictionary with the whitelist
        try:
          with open(self.whitelist, 'r') as f:
              for bc_sequence in f:
                  bc_sequence = bc_sequence.strip()
                  if len(bc_sequence) != self.bc_length:
                      raise ValueError("ERROR[barcodes_10xG]\tWhitelist barcode length" +
                                       " does not match barcode trim length")
                  if bc_sequence in self._gbcDict:
                      raise ValueError("ERROR[barcodes_10xG]\tDuplicate barcode" +
                                       " found in the whitelist")
                  self._gbcDict.add(bc_sequence)
        except EnvironmentError:
            sys.stderr.write('ERROR[barcodes_10xG]\tError reading in barcode whitelist')
            raise
        self.wl_size = len(self._gbcDict)
        if self.verbose:
            sys.stderr.write(("PROCESS\tNOTE\tFinished reading in {0} barcodes" +
                             " from whitelist.\n").format(self.wl_size))

    def __checkHammingOne(self, sequence):
        res = []
        i = 0
        try:
            while i < len(sequence):
                for j in self.__basedict.get(sequence[i]):
                    testSeq = sequence[:i] + j + sequence[i + 1:]
                    if testSeq in self._gbcDict:
                        res.append(testSeq)
                i += 1
        except KeyError:
            sys.stderr.write("ERROR[barcodes_10xG]\tUnknown base character sequence found")
            raise
        return res

    def check_barcode(self, barcode):
        """
        for an extracted barcode, compare it to the whitelist and
        return whitelisted barcode and status.
        """
        status = "UNKNOWN"
        if 'N' not in barcode and barcode in self._gbcDict:  # barcode matches whitelist
            status = "MATCH"
            gem_bc = barcode
        else:
            hamming = self.__checkHammingOne(barcode)
            if len(hamming) == 0:  # greater than 1 hamming distance
                status = "UNKNOWN"
                gem_bc = None
            elif len(hamming) == 1:  # single hit hamming distance of 1
                status = "MISMATCH1"
                gem_bc = hamming[0]
            else:  # multihit hamming distance of 1
                status = "AMBIGUOUS"
                gem_bc = None
        self._gbcCounter[(barcode, gem_bc)] += 1
        self.statusCounter[status] += 1
        return gem_bc, status

    def get_gbcCounter_size(self):
        return len(self._gbcCounter)

    def get_gbcCounter_items(self):
        return self._gbcCounter.values()

    def print_gbcCounter(self, filename):
        """
        Function to print to file the barcode counter information
        """
        try:
            with open(filename, 'w') as f:
                [f.write('{0}\t{1}\t{2}\n'.format(key[0], key[1], value)) for key, value in self._gbcCounter.items()]
         except EnvironmentError:
            sys.stderr.write('ERROR[barcodes_10xG]\tError writing in barcode counter')
            raise
