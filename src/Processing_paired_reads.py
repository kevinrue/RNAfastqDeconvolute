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
        (Format to be defined, depending on how we decide to scan the reads for the adapter sequences. For instance,hep
        Nick required a STRING "I3S3D0" because he used the Perl module StringApprox which requires exactly this format
        to do pattern matching. All Nick had to do is forward the STRING to the StringApprox module.)

    -I, -illumina_version FLOAT
        Version of the Illumina Genome Analyser pipeline used to generate the fastq files. This will affect the
        conversion of the quality score from ASCII to PHRED using an offset of 33 or 64.

    -p, --percent_max FLOAT
        Maximal percentage of the read length allowed below the specified PHRED threshold (see -P, --phred). Reads with
        a larger observed percentage will be filtered out due to "poor overall quality".

    -P, --phred INTEGER
        Threshold used to define sequenced nucleotides as "poor quality" (see -p, --percent_max).


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
# Module argparse allows to process the arguments from the command line automatically instead of manually parsing the
# individual elements.
# argparse will return an error if a required argment is missing.
# argparse will convert the value provided in the command line to the data type specified in the parser definition.
import argparse
# Module os.path allows to check if paths and files exist.
import os.path
# Custom module to read compressed fastq files. Read it as "from folder RNAseqIO, import the script Parser".
# The RNAseqIO folder needs to be on the $PYTHONPATH or in the folder of the script initially called from the command
# line.
from RNAseqIO import Parser


def main():
    # Define the argument parser
    parser = argparse.ArgumentParser(description="This package deconvolutes PE Illumina reads")
    # List the mandatory options in a separate section "essential arguments" of the help message
    essential_options = parser.add_argument_group('essential arguments')
    # Example of a required argument to save as STRING variable (default)
    # metavar affects the usage help message
    # dest is the name of the local variable where the filename will be saved "args.<dest>"
    essential_options.add_argument('-f', '--forward', required=True,
                                   help='GZIP-compressed fastq file containing the forward reads.',
                                   metavar='pe1.fastq.gz',
                                   dest='forward_file')
    # Example of a required argument to save as FLOAT number
    essential_options.add_argument('-I', '--illumina_version', required=True, type=float,
                                   help='Version of the Illumina Genome Analyser pipeline used to generate the fastq \
                                   files. This will affect the conversion of the quality score from ASCII to PHRED \
                                   score using an offset of 33 or 64.', metavar='x.x')
    # Lists the optional arguments in the default section named "optional arguments:"
    # Example of an optional argument which will be set to 25 if not specified, and saved to args.percent_max
    parser.add_argument('-p', '--percent_max', type=float, default=25,
                        help="Maximal percentage of trimmed nucleotides allowed below the Phred threshold ()."
                             "(see -P, --phred). Default is 25.", metavar='p')
    # Example of an optional argument which will be set to 20 if not specified and saved to args.phred
    parser.add_argument('-P', '--phred', type=int, default=20,
                        help="Minimal accepted Phred score used to define sequenced nucleotides as \"good quality\""
                             "(see -p, --percent_max). Default is 20.", metavar='P')
    # parse command line options according to the rules defined above
    args = parser.parse_args(sys.argv[1:])
    # The elements in the agrs variable will be printed beside Namespace once the code is run
    # For training purpose: print all the arguments found TODO: remove in final code
    print(args)
    # For training purpose: print here  how to access the value of the option "forward" TODO: remove in final code
    #print args.forward_file
    # Check that the forward_file provided does exist
    if not os.path.isfile(args.forward_file):
        print("Error: File of forward reads was not found: %s" % args.forward_file)
        sys.exit(2)
    # Opens the forward read file
    forward_parser = Parser.FastqgzParser(args.forward_file)
    # Get the next read
    read = forward_parser.nextRead()
    # Initialises a counter storing the number of read (pairs) filtered out because of adapter contamination
    adapter_count = 0
    # Converts the user-defined Phred threshold to a Ascii-compatible value
    if 1 <= args.illumina_version < 1.8:
        ascii_phred_threshold = args.phred + 64 - 1
    elif args.illumina_version >= 1.8:
        ascii_phred_threshold = args.phred + 33 - 1
    else:
        print("This illumina version: %d is not supported; please check again your illumina version for phred score \
               encoding or seek advice about this script!\n" % args.illumina_version)
        sys.exit(3)
    # For training purpose: print the value TODO: remove in final code
    print(
        "ascii_phred_threshold: %i (Phred: %i + Illumina v%s offset: %i - 1 to avoid doing -1 for each read)" %
        (ascii_phred_threshold, args.phred, args.illumina_version, ascii_phred_threshold - args.phred + 1))
    # While the last read parsed is not empty (= end of file not reached), process the read
    while read.header_line:
        print(read)
        # TODO replace the numbers in the line below by numbers calculated from the parsed command line
        read.trim(1, 89)  # trim the first and last bases
        print(read)
        # 20 the Phred threshold for testing here, 64 the offset for Illumina 1.5, and -1 for mathematical reasons (the
        # define_quality_status function uses percentile to check the quality much faster than a per-base counter)
        read.define_quality_status(ascii_phred_threshold, args.percent_max)
        # For training purpose: print quality_status attribute (True if accepted quality) TODO: remove in final code
        print("read.quality_status: %s" % read.quality_status)
        # Check whether the adapter sequence is present without mismatch
        # TODO replace the values in the line below by values parsed from the adaptor file and command line
        read.define_adapter_presence("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", 3)
        # For training purpose: print quality_status attribute (True if adapter present) TODO: remove in final code
        print("read.adapter_present: %s" % read.adapter_present)
        # adds one to adapter counter in adapter is present
        if read.adapter_present:
            adapter_count += 1
        # Moves on to the next (we don't want to process the same read eternally, do we?)
        read = forward_parser.nextRead()
    # For training purpose: print the count of reads with adapter detected TODO: remove in final code
    print("adapter_count: %i" % adapter_count)
    # Close the file stream
    forward_parser.close()


if __name__ == "__main__":
    main()
