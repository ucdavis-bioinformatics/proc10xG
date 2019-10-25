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

def sp_gzip_read(file,bufsize):
    p = subprocess.Popen(["gzip", "--decompress", "--to-stdout", file],
                         stdout=subprocess.PIPE,
                         stderr=None,
                         bufsize=bufsize,
                         text=True,
                         preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    if p.returncode:
        raise
    return p.stdout


def sp_gzip_write(file, bufsize=-1):
    filep = open(file, 'w')
    p = subprocess.Popen('gzip',
                         stdin=subprocess.PIPE,
                         stdout=filep,
                         shell=True,
                         bufsize=bufsize)
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
