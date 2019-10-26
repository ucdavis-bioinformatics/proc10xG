import os
import errno
import string
import numpy
import subprocess

def median(lst):
    return numpy.median(numpy.array(lst))

'''
Gzip utilities, run gzip in a subprocess
'''
# unversal_newlines is being used here for backwards compatability to
# python 3.6
def sp_gzip_read(file,bufsize=-1):
    p = subprocess.Popen(["gzip", "--decompress", "--to-stdout", file],
                         stdout=subprocess.PIPE,
                         stderr=None,
                         bufsize=bufsize,
                         universal_newlines=True)
    if p.returncode:
        raise
    return p.stdout


# unversal_newlines is being used here for backwards compatability to
# python 3.6
def sp_gzip_write(file, bufsize=-1):
    filep = open(file, 'w')
    p = subprocess.Popen('gzip',
                         stdin=subprocess.PIPE,
                         stdout=filep,
                         shell=True,
                         bufsize=bufsize,
                         universal_newlines=True)
    return p.stdin


def make_sure_path_exists(path):
    """
    Try and create a path, if not error
    """
    if path != '':
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    return path


rcs = str.maketrans('TAGCtagc', 'ATCGATCG')


def revcomp(seq):
    return seq.translate(rcs)[::-1]


def rev(seq):
    return seq[::-1]
