#!/usr/bin/env python

"""
NAME
    Processing_paired_reads.py - Deconvolute, trims and filters fastq files

SYNOPSIS
    Processing_paired_reads.py --forward raw1.fastq.gz --reverse raw2.fastq.gz --index indices.txt ...

DESCRIPTION
    Parses two fastq files containing paired reads from a sequencing lane, deconvolutes the reads according to
    a barcode/sample mapping file, filters out reads contaminated with adapter sequence provided by the user, filters
    out reads of poor quality based on user-defined thresholds, and trims the reads according to user instructions.

    -f, --forward file.fastq.gz
        GZIP-compressed fastq file containing the forward reads.

    -r, --reverse file.fastq.gz
        GZIP-compressed fastq file containing the reverse reads.

    -i, --index file.txt
        TAB-separated file containing two columns. First, the barcode as expected in the read; Second, an identifier
        unique to the corresponding sample (e.g. animal_time_treatment).

    -t, --trim STRING
        The trimming parameters, specific in the format "5end_INTEGER1_3end_INTEGER2". INTEGER1 bases will be trimmed
        off the 5' end and vice versa.

    -a, --adapter file.txt
        TAB-separated file containing three columns. First, the adapter sequence as expected in the read; Second, an
        identifier unique to the adapter sequence; Third, one of "forward" or "reverse" specifying the mate in which
        to look for the adapter sequence.

    -A, --adapter_param STRING
        (Format to be defined, depending on how we decide to scan the reads for the adapter sequences. For instance,
        Nick required a STRING "I3S3D0" because he used the Perl module StringApprox which requires exactly this format
        to do pattern matching. All Nick had to do is forward the STRING to the StringApprox module.)

    -I, -illumina_version DECIMAL
        Version of the Illumina Genome Analyser pipeline used to generate the fastq files. This will affect the
        conversion of the quality score from ASCII to PHRED using an offset of 33 or 64.

    -p, --percent _thresh FLOAT
        Maximal percentage of the read length allowed below the specified PHRED threshold (see -P, --phred). Reads with
        a larger observed percentage will be filtered out due to "poor overall quality".

    -P, --phred INTEGER
        Threshold used to define sequenced nucleotides as "poor quality" (see -p, percent_thresh).


The main function was based on a template posted at https://www.artima.com/weblogs/viewpost.jsp?thread=4829
by Guido van Rossum (creator of Python).

The overall deconvolution procedure programmed here was based on the original Perl version programmed by Nicolas
Nalpas and Stephen Park (Animal Genomics Lab, UCD, Dublin, Ireland).
"""

__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

# Module sys allows to interact with the operating system
# for instance, collect flags and arguments from the command line
import sys
# Module getopt  allows to process the arguments from the command line
# automatically instead of manually parsing the individual elements
import getopt


def main():
    # parse command line options
    try:
        # saves the expected option/argument pairs in opts, and the remaining arguments in args
        opts, args = getopt.getopt(sys.argv[1:], 'hf:r:i:t:a:A:I:p:P:',
                                   ["help", "forward=", "reverse=", "index=", "trim=", "adapter=",
                                    "adapter_param=", "illumina_version=", "percent_thresh=",
                                    "phred="])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        # For training purpose: shows all the options identified
        print(o, a)
        # if the "help" option is given, prints the docstring (above)
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
    # process arguments
    # arguments are what is left after all the expected option have been parsed
    for arg in args:
        print(arg)  # process() must be defined elsewhere


if __name__ == "__main__":
    main()
