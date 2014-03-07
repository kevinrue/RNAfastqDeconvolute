#!/usr/bin/env python

"""
NAME
    Processing_paired_reads.py - Deconvolute, trims and filters fastq files

SNYOPSIS
    Processing_paired_reads --forward raw1.fastq.gz --reverse raw2.fastq.gz --index indices.txt ...

DESCRIPTION
    Parses two fastq files containing paired reads from a sequencing lane, deconvolutes the reads according to
    a barcode/sample mapping file, filters out reads contaminated with adapter sequence provided by the user, filters
    out reads of poor quality based on user-defined thresholds, and trims the reads according to user instructions.

    -a, --adapter file
        TAB-separated file containing three columns. First, the adapter sequence as expected in the read; Second, an
        identifier unique to the adapter sequence; Third, one of "forward" or "reverse" specifying the mate in which
        to look for the adapter sequence.

This is Contained in variable named __doc__

This code was based on an example found at https://www.artima.com/weblogs/viewpost.jsp?thread=4829
by Guido van Rossum (creator of Python).
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
        opts, args = getopt.getopt(sys.argv[1:], 'hf:r:i:I:t:a:A:p:P:',
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
