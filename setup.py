#!/usr/bin/env python3

from setuptools import setup

config = \
    {
        'name': 'proc10xG',
        'packages': ['proc10xG'],
        'version': '0.1.0',
        'description': 'Processing of 10X Genomics projects',
        'author': 'Matt Settles',
        'author_email': 'settles@ucdavis.edu',
        'url': 'https://github.com/ucdavis-bioinformatics/proc10xG',
        'download_url': 'https://github.com/ucdavis-bioinformatics/proc10xG',
        'keywords': ["Bioinformatics", "DNA Sequencing", "10X Genomics"],
        'classifiers': [
            "Programming Language :: Python",
            "Programming Language :: Python :: 3.6",
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: Apache Software License",
            "Operating System :: OS Independent",
            "Topic :: Software Development :: Libraries :: Python Modules",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
            ],
        'install_requires': ['numpy'],
        'entry_points': {
          'console_scripts': [
              'process_10xReads = proc10xG.process_10xReads:main']
              },
        'data_files': [('data/barcodes', ['data/barcodes/4M-with-alts-february-2016.txt',
                                          'data/barcodes/737K-april-2014_rc.txt',
                                          'data/barcodes/737K-august-2016.txt'])]
    }

setup(**config)
