#!/usr/bin/env python

"""docstring:
This serves as a long usage message.

This script is the main script for the processing of paired reads in fastq

This code was based on an example found at https://www.artima.com/weblogs/viewpost.jsp?thread=4829
by Guido van Rossum (creator of Python).
"""

__author__ = 'David Magee, Carolina Correia, and Kevin Rue'
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
        opts, args = getopt.getopt(sys.argv[1:], 'ha:b:c:d:', ["help"])
    except getopt.error, msg:
        print msg
        print "for help use --help"
        sys.exit(2)
    # process options
    for o, a in opts:
        # if the "help" option is given, prints the docstring (above)
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)
        # For training purpose: shows all the options identified
        print(o, a)
    # process arguments
    # arguments are what is left after all the expected option have been parsed
    for arg in args:
        print(arg)  # process() must be defined elsewhere


if __name__ == "__main__":
    main()
