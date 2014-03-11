__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module os allows to access OS-specific line separator.
import os
# Module scipy allows to calculate quantiles, useful to define the PHRED score below which a given percentage of the
# bases in a read are found.
import scipy
# Module fuzzysearch allows to find approximate matches, useful to detect adaptor sequences in reads with possible
# sequencing errors.
import fuzzysearch


class Read:
    """Sequenced read and associated data"""

    def __init__(self, header_line, sequence_line, separator_line, quality_line):
        """Constructor for Read"""
        self.header_line = header_line
        self.sequence_line = sequence_line
        self.separator_line = separator_line
        self.quality_line = quality_line
        # If acceptable quality then set to True, otherwise set to False
        # Defaults to True (if no quality check is done, quality is assumed to be good)
        self.quality_status = True
        # If adapter is present then set to True, otherwise set to False
        # Defaults to False (if no adapter check is done, it is assumed to be absent)
        self.adapter_present = True

    def trim(self, left, right):
        """Crops the read sequence and associated quality_line keeping positions five_start to three_end.
        Be careful to give the positions in Python position numbering (starts at 0). The pain of calculating the proper
        coordinates before calling the function will save a lot of time avoiding re-calculating them every time the
        function is called.

        :rtype : None
        Args:
            self, left, right

        Returns:
            None
        """
        self.sequence_line = self.sequence_line[left:right]
        self.quality_line = self.quality_line[left:right]

    def define_quality_status(self, threshold, percentage):
        """Sets the quality_status attribute according to whether stricly more than percentage of bases in the read are
         strictly below phred.
         Be careful to give a threshold that is the desired minimal Phred score -1. The reason can be demonstrated with
         the extreme example where the read has the minimal number of bases at the minimal score or higher. The pain of
         calculating the proper coordinates before calling the function will save a lot of time avoiding re-calculating
         them every time the function is called.

        :rtype : None
        Args:
            self, phred, percentage

        Returns:
            None
        """
        # Converts all the ascii caracters in the quality line into their ascii code
        quality_ascii = [ord(c) for c in self.quality_line]
        # Defines whether the read is acceptable (True) or not (False)
        self.quality_status = scipy.percentile(quality_ascii, percentage) > threshold
        # For training purpose: prints the value of the percentile used for threshold TODO remove in final code
        print scipy.percentile(quality_ascii, percentage)


    def define_adapter_presence(self, adapter, edit_threshold):
        """Set the adapter_absent attribute according to whether a match is found with an edit distance less or equal to
        edit_threshold.

        :rtype : None
        Args:
            self, adapter, edit_threshold

        Returns:
            None
        """
        self.adapter_present = len(
            fuzzysearch.find_near_matches_with_ngrams(adapter, self.sequence_line, edit_threshold)) > 0


    def __str__(self):
        """This function defines the string representation displayed when calling the print function on a Read object"""
        # Uses the nw-line separator of the operating system running this script to separate the different elements of
        # the read on different lines at the screen.
        return os.linesep.join([self.header_line, self.sequence_line, self.separator_line, self.quality_line])

