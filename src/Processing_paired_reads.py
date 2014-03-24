#!/usr/bin/env python

"""
COMMENTS ON RNA-SEQ READS FROM BGI
    The header contains the index (the first six nucleotides of the read). Some indexes contain eight nucleotides with
    the last two nucleotides added by the BGI. Only use the first six nucleotides for the indexed matching of the
    libraries.

    The reverse reads are NOT reverse complemented by the BGI. However, the index of the forward and reverse reads
    from the same mate pair have the same index sequence (this allows mate pairs to be identified)

    If adapter detected in forward reads, then adapter must be detected in the reverse reads. The forward adapter
    sequence is the same as the Illuumina TruSeq_Adapter (AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC) and the reverse adapter
    is Illumina_Single_End_PCR_Primer_1 (AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA). The code will need to be amended to
    take care of mismatches within the adapter.

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

    -b, --barcodes file.txt
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

# Module argparse allows to process the arguments from the command line automatically instead of manually parsing the
# individual elements.
# argparse will return an error if a required argment is missing.
# argparse will convert the value provided in the command line to the data type specified in the parser definition.
import argparse
# Module collections allows to store multiple variable in a named tuple. Very clean for start and stop trimming indices
import collections
# Module datetime allows to print the current time at the start and the end of the script
import datetime
# Module os allows to dynamically use the OS-dependent newline character
import os
# Module os.path allows to check if paths and files exist. It is also useful to print the current working directory.
import os.path
# Module sys allows to interact with the operating system
# for instance, collect flags and arguments from the command line
import sys
# Custom module to log read deconvolution stastistics and write them to a report file
import Loggers
# Custom module to read compressed fastq files. Read it as "from folder RNAseqIO, import the script Parser".
# The RNAseqIO folder needs to be on the $PYTHONPATH or in the folder of the script initially called from the command
# line.
# Parser is a module generated by Kevin (from the Parsers.py script) that defines Fastqgz Parser that is called below.
# Make sure that the RNAseqIO folder path is set correctly. To do this, right click on the src folder at the top
# right hand side of the page and select 'Mark Directory As'. Then select 'as source root). This resolved the
# previous error message stating that RNASeqIO was an unresolved reference.
import Parsers


def main():
    # Informative message: print the starting time (strip off the milliseconds)
    print("Info: Start time:", format(datetime.datetime.now().replace(microsecond=0)))
    # Informative message: print the current path
    print("Info: Working directory: {0:s}".format(os.getcwd()))

    # Define the argument parser
    # Please follow an alphabetical order for ease of code searching
    parser = argparse.ArgumentParser(description="This package can deconvolute and filter PE Illumina reads.")
    # List the mandatory options in a separate section "essential arguments" of the help message
    essential_options = parser.add_argument_group('mandatory arguments')
    # Example of a required argument to save as STRING variable (default)
    # metavar affects the usage help message
    # dest is the name of the local variable where the filename will be saved "args.<dest>"
    essential_options.add_argument('-f', '--forward', required=True,
                                   help='GZIP-compressed fastq file containing the forward reads.',
                                   metavar='pe1.fastq.gz',
                                   dest='forward_file')
    # Example of a required argument to save as FLOAT number
    essential_options.add_argument('-r', '--reverse', required=True,
                                   help='GZIP-compressed fastq file containing the reverse reads.',
                                   metavar='pe2.fastq.gz',
                                   dest='reverse_file')
    # Required argument to provide barcodes for deconvolution
    essential_options.add_argument('-b', '--barcodes',
                                   help="TAB-separated file mapping barcodes sequences (1st column) to samples (2nd "
                                        "column)."
                                        "No file skips the deconvolution.", metavar='barcodes.txt')
    # Example of a required argument to save as FLOAT number
    essential_options.add_argument('-I', '--illumina_version', required=True, type=float,
                                   help='Version of the Illumina Genome Analyser pipeline used to generate the fastq \
                                   files. This will affect the conversion of the quality score from ASCII to PHRED \
                                   score using an appropriate offset of 33 or 64.', metavar='x.x')

    # Lists the optional arguments in the default section named "optional arguments:"
    # Optional argument to provide adapters for filtering
    parser.add_argument('-a', '--adapters', default=None,
                        help="TAB-separated file mapping adapter sequences (1st column) to unique identifiers (2nd "
                             "column), and mate {forward, reverse} to screen for the adapter (3rd column). No file "
                             "skips"
                             " the filtering.", metavar='adapters.txt')
    # Example of an optional argument which will be set to 25 if not specified, and saved to args.percent_max
    parser.add_argument('-p', '--percent_max', type=float, default=25,
                        help="Maximal percentage of trimmed nucleotides allowed below the Phred threshold "
                             "(see -P, --phred). Default is 25.", metavar='p')
    # Example of an optional argument which will be set to 20 if not specified and saved to args.phred
    parser.add_argument('-P', '--phred', type=int, default=None,
                        help="Minimal accepted Phred score used to define sequenced nucleotides as \"good quality\" "
                             "(see -p, --percent_max). Omitting the option skips quality filtering.", metavar='P')
    # Optional argument to trim the read, 2 values expected after the flag for 5' and 3' trimmed length,
    parser.add_argument('-t', '--trim', default=None, nargs=2, type=int,
                        help='Number of bases to trim from the 5\' and 3\' of the read sequence and associated quality '
                             'string. Example "-t 4 3" will trim the first 4 bases and the last 3 bases. Omitting the '
                             'option skips the trimming.', metavar="T")
    # parse command line options according to the rules defined above
    args = parser.parse_args(sys.argv[1:])
    # The elements in the agrs variable will be printed beside Namespace once the code is run
    # For training purpose: print all the arguments found TODO: remove in final code
    #print("Test: args: %s\n" % args)

    # Barcode is a required argument, if we reach here, it means the user gave a value
    # Informative message
    print("Info: Barcode deconvolution: using file %s" % args.barcodes)
    # If the file does not exist, print an error message and exit the script
    if not os.path.isfile(args.barcodes):
        print("Error: File of barcodes was not found: %s" % args.forward_file)
        sys.exit(4)
    # If we reach this line onward, the file does exist
    # Create a BarcodesParser object
    barcode_parser = Parsers.BarcodesParser(args.barcodes)
    # Make the BarcodesParser import the mapping of expected barcodes and corresponding samples
    # No need to store the output in a variable. The BarcodesParser does not return anything, it stores the
    # information in one of its own variable accessible here as "barcode_parser.expected"
    barcode_parser.parse_barcode_file()
    # For training purposes, print the dictionary mapping TODO: remove in final code
    #print("Test: barcode_parser.expected: %s" % barcode_parser.expected)
    # From the first barcode, define the expected barcode length
    barcode_length = len(list(barcode_parser.expected.keys())[0])
    # Informative message
    print("Info: Expected barcode length: %i (observed from a sample in barcode file)" % barcode_length)

    # If the user gave an adapter file
    if args.adapters:
        # Informative message
        print("Info: Adapter filtering: using file %s" % args.adapters)
        # If the file does not exist, print an error message and exit the script
        if not os.path.isfile(args.adapters):
            print("Error: File of adapters was not found: %s" % args.adapters)
            sys.exit(4)
        # If we reach this line onward, the file does exist
        # Create an AdapterParser object
        adapter_parser = Parsers.AdapterParser(args.adapters)
        # Make the AdapterParser import the forward and reverse adapters and associated information TODO remove in final
        # code
        adapters = adapter_parser.parse_adapters_file()
        #print("Test: adapters:", adapters, '\n')
        print("Info: Adapter searched in forward reads: %s" % adapters.forward)
        print("Info: Adapter searched in reverse reads: %s" % adapters.reverse)
    # if the user did not provide a barcode indexing file
    else:
        # Informative message
        print("Info: Adapter filtering: skipped")
        #set the adapters as None
        adapters = None

    # If the user gave a phred threshold
    if args.phred:
        # Informative message
        print(
            'Info: Quality filtering: Maximum of {0:.3f}% bases in each mate are allowed strictly below {1:d} Phred '
            'score.'.format(args.percent_max, args.phred))
        # Converts the user-defined Phred threshold to a Ascii-compatible value
        if 1 <= args.illumina_version < 1.8:
            ascii_phred_threshold = args.phred + 64 - 1
        elif args.illumina_version >= 1.8:
            ascii_phred_threshold = args.phred + 33 - 1
        else:
            print(
                "Error: The Illumina version: %.1f is not supported; please check again your illumina version for "
                "phred score encoding or seek advice about this script!%s" % (args.illumina_version, os.linesep))
            sys.exit(4)
            # For training purpose: print the value TODO: remove in final code
            #print(
            #    "Test: ascii_phred_threshold: %i (Phred: %i + Illumina v%s offset: %i - 1 to avoid doing -1 for each
            #  read "
            #    "in the testing step)\n" % (
            #        ascii_phred_threshold, args.phred, args.illumina_version, ascii_phred_threshold - args.phred + 1))
    # if the user did not provide a phred threshold
    else:
        # Informative message
        print("Info: Quality filtering: skipped")
        # leave the args.phred value as None to skip quality filtering
        # set ascii_phred_threshold to avoid a PyCharm warning message about possibly undefined variable later
        ascii_phred_threshold = None

    # Check that the forward_file provided does exist
    if not os.path.isfile(args.forward_file):
        print("Error: File of forward reads was not found: %s" % args.forward_file)
        sys.exit(2)
    # Opens the forward read file using the Parser function (generated from the Parsers.py script)
    forward_parser = Parsers.FastqgzParser(args.forward_file)
    # Get the next read (first) as defined in the Parser function (generated by Parsers.py). The Parser function
    # calls each read as a series of four lines (one line for the header, sequence, separator and quality)
    forward_read = forward_parser.next_read()
    # check that the reverse file provided does exist
    if not os.path.isfile(args.reverse_file):
        print("Error: File of reverse reads was not found: %s" % args.reverse_file)
        sys.exit(3)
    # Opens the reverse read file using the Parser function (generated from the Parsers.py script)
    reverse_parser = Parsers.FastqgzParser(args.reverse_file)
    # Get the next read (first) as defined in the Parser function (generated by Parsers.py). The Parser function
    # calls each read as a series of four lines (one line for the header, sequence, separator and quality)
    reverse_read = reverse_parser.next_read()
    # The Parser can be run on both the forward_read and reverse_read variables

    # Now that we have the first read, we can do a few more checks assuming that the first read is representative
    # of all the reads. This assumption is reasonable regarding read length, header structure, ..
    # Return an error and exit if the first read length is 0
    if len(forward_read.sequence_line) == 0:
        print("Error: Read length seems to be 0 (observed from the first read)")
    # Informative message:
    print("Info: Read length seems to be %i (observed from the first read)" % len(forward_read.sequence_line))
    # If the user gave trimming values
    # Check that the trimmed read will have at least 1 base left
    if args.trim:
        # Return an error and exit if the sum of trimming length is larger than the read length
        if sum(args.trim) >= len(forward_read.sequence_line):
            print(
                "Error: The trimming parameters add up to equal or more bases than there are in the first read. Please "
                "check your data or seek advice about this script!")
            sys.exit(4)
        # If we reach here, we know that the trimming will leave 1 or more bases in the read
        # if the user explicitly asked to trim nothing at both end
        if args.trim == [0, 0]:
            # set args.trim to None, this will skip trimming
            args.trim = None
            # Informative message
            print("Info: Trimming: skipped (0, 0)")
        # if the user really asked to trim something at either end
        else:
            # calculate the Python start and stop index once and for all (saves the time of doing it for each read)
            # Prepare a named tuple structure to cleanly return the two adapters in one variable
            trim = collections.namedtuple('TrimIndices', ['start', 'stop'])
            # Fill in the tuple structure with the adapter objects and return it
            trim_indices = trim(start=args.trim[0], stop=len(forward_read.sequence_line) - args.trim[1])
            # Informative message
            print(
                "Info: Trimming: %i bases will be trimmed from the left and %i bases will be trimmed from the right, "
                "leaving %i bases in the trimmed read." % (
                    args.trim[0], args.trim[0], len(forward_read.sequence_line) - sum(args.trim)))
    # if the user did not provide a phred threshold
    else:  # Informative message
        print("Info: Trimming: skipped")
    # leave the args.trim value as None

    # Initialises a ReadLogger to store the deconvolution statistics
    read_logger = Loggers.ReadLogger(args.forward_file, barcode_parser)

    # Parses the reads and process them
    # While the last read parsed is not empty (= end of file not reached), process the read
    # We will assume that if the header line is not empty, the other fields are not either
    # Also we only check the first read to save more time (a single test per mate pair)
    while forward_read.header_line:
        # For training purpose print the
        #print("Test: forward_read: %s\n" % forward_read)
        #print("Test: reverse_read: %s\n" % reverse_read)
        # run the while function on both the forward and reverse read simultaneously

        # DECONVOLUTION #
        # Only perform deconvolution step if barcodes were given
        if args.barcodes:
            # Define which sample the raed should be assigned to
            barcode_parser.assign_read_to_sample(forward_read, barcode_length)
            #print("Test: read.sample: %s\n" % forward_read.sample)
            # if the barcode was assigned to any sample
            if forward_read.sample:
                # log it
                read_logger.assigned[forward_read.sample] += 1
            # if the barcode was not assigned to any sample
            else:
                # log it
                read_logger.unassigned += 1
                # and no need to further process it, move to the next read
                forward_read = forward_parser.next_read()
                reverse_read = reverse_parser.next_read()
                # and skip the rest of this loop (which is why we had to read the next read here)
                continue

        # TRIMMING #
        # Only perform the trimming step if trimming values were given and at least one of them is different from 0
        if args.trim:
            forward_read.trim(trim_indices.start, trim_indices.stop)
            # Trim is a function defined in the SeqDataTypes.py script
            # In pycharm, to see what Trim does hover over the word 'trim' while hitting CTRL
            reverse_read.trim(trim_indices.start, trim_indices.stop)  # trim the reverse read
            #print("Test: forward_read after trimming: %s\n" % forward_read)
            #print("Test: reverse_read after trimming: %s\n" % reverse_read)

        # ADAPTER #
        # Only perform the adapter detection step if the adapter was given at the command line
        if args.adapters:
            # TODO replace the 3 values in the line below by values parsed from the adaptor file and command
            # line
            forward_read.define_adapter_presence_substitutions_only(adapters.forward.sequence, 3)
            # if the adapter was found in the forward read
            if forward_read.adapter_present:
                # log it
                read_logger.adapter_excluded[forward_read.sample] += 1
                # and no need to further process it, move to the next read
                forward_read = forward_parser.next_read()
                reverse_read = reverse_parser.next_read()
                # and skip the rest of this loop (which is why we had to read the next read here)
                continue
            # if we reach here the adapter was not found in the forward mate, check the reverse mate
            reverse_read.define_adapter_presence_substitutions_only(adapters.reverse.sequence, 3)
            # if the adapter was found in the forward read
            if reverse_read.adapter_present:
                # log it
                read_logger.adapter_excluded[forward_read.sample] += 1
                # and no need to further process it, move to the next read
                forward_read = forward_parser.next_read()
                reverse_read = reverse_parser.next_read()
                # and skip the rest of this loop (which is why we had to read the next read here)
                continue
                # if we reach here the adapter was absent from both mates
                # For training purpose: print quality_status attribute (True if adapter present) TODO: remove in
                # final code
                #print("Test: forward_read.adapter_present: %s\n" % forward_read.adapter_present)
                #print("Test: reverse_read.adapter_present: %s\n" % reverse_read.adapter_present)

        # QUALITY #
        # Only perform the quality filtering step if a Phred threshold was provided (not None, see above)
        if args.phred:
            # define_quality_status function uses percentile to check the quality much faster than a per-base
            #  counter)
            forward_read.define_quality_status(ascii_phred_threshold, args.percent_max)
            # if the forward read is poor quality
            if not forward_read.quality_status:
                # log it
                read_logger.quality_excluded[forward_read.sample] += 1
                # and no need to further process it, move to the next read
                forward_read = forward_parser.next_read()
                reverse_read = reverse_parser.next_read()
                # and skip the rest of this loop (which is why we had to read the next read here)
                continue
            # if we reach here the quality of the forward mate was acceptable, check the reverse mate
            reverse_read.define_quality_status(ascii_phred_threshold, args.percent_max)
            # if the forward read is poor quality
            if not reverse_read.quality_status:
                # log it
                read_logger.quality_excluded[forward_read.sample] += 1
                # and no need to further process it, move to the next read
                forward_read = forward_parser.next_read()
                reverse_read = reverse_parser.next_read()
                # and skip the rest of this loop (which is why we had to read the next read here)
                continue
                # if we reach here the quality of both mates was acceptable
                # For training purpose: print quality_status attribute (True if accepted quality) TODO: remove in
                #  final code
                #print("Test: read.quality_status: %s\n" % forward_read.quality_status)

        # if we reach here, the read pair succesfully passed deconvolution and filtering
        # log it
        read_logger.accepted[forward_read.sample] += 1

        # Moves on to the next (we don't want to process the same read eternally, do we?)
        forward_read = forward_parser.next_read()
        reverse_read = reverse_parser.next_read()

    # For training purpose: print the stats recorded in the read.logger TODO: remove in final code
    #print("Test: read_logger.assigned: ", read_logger.assigned)
    #print("Test: read_logger.unassigned: %i" % read_logger.unassigned)
    #print("Test: read_logger.quality_excluded:", read_logger.quality_excluded)
    #print("Test: read_logger.adapter_excluded:", read_logger.adapter_excluded)
    print('Info: Total reads processed: {0:,}'.format(sum(read_logger.assigned.values()) + read_logger.unassigned))

    # Close the file stream of the raw read files
    forward_parser.close()
    reverse_parser.close()

    # Write the deconvolution statistics in the report file
    read_logger.write_stats()

    # Informative message
    print("Info: End time:", datetime.datetime.now().replace(microsecond=0))


if __name__ == "__main__":
    main()
