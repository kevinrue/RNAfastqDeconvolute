__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module scipy allows to calculate quantiles, useful to define the PHRED score below which a given percentage of the
# bases in a read are found.
import scipy
# Custom Module which contains a few functions for approximate matching
import ApproxMatch


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
        self.adapter_present = False

    def trim(self, start, stop):
        """Crops the read sequence and associated quality_line keeping positions start to stop.
        Be careful to give the positions in Python position numbering (starts at 0). The pain of calculating the proper
        coordinates before calling the function will save a lot of time avoiding re-calculating them every time the
        function is called.

        Args:
            self, start, stop

        Returns:
            None
        """
        # Trims the sequence line to the desired window
        self.sequence_line = self.sequence_line[start:stop]
        # Trims the quality line to the same window
        self.quality_line = self.quality_line[start:stop]

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
        #print('Test: percentile ({0}%): {1:.3f}\n'.format(percentage, scipy.percentile(quality_ascii, percentage)))

    def define_adapter_presence_substitutions_only(self, adapter, max_substitutions):
        """Sets the adapter_absent attribute according to whether a match is found with a number of substitutions
        less or equal to
        edit_threshold.

        :rtype : None
        Args:
            self, adapter, edit_threshold

        Returns:
            None

        Args:
            self, adapter, max_substitutions

        Returns:
            None
        """
        # for each subsequence of the sequence_line of length identical to the adapter (last starting position is
        # length of the adapter before the end of the read)
        for index in range(len(self.sequence_line) - len(adapter) + 1):
            # check if the adapter is less than max_substitutions away from the subsequence
            result = ApproxMatch.approx_substitute(adapter, self.sequence_line[index:index + len(adapter)],
                                                   max_substitutions)
            # If a match was found (result is TRUE)
            if result:
                # set the adapter presence to True (and stop the function here)
                self.adapter_present = True
                return None
                # If no substring of sequence_line approximately matches the adapter, the function will eventually
                # leave the
                # loop above.
                # The simple fact of arriving here proves that no match was found, therefore leave the adapter
                # preence to False
                # and simply leave the function

    def __str__(self):
        """This function defines the string representation displayed when calling the print function on a Read object"""
        # Uses the new-line separator of the operating system running this script to separate the different elements of
        # the read on different lines at the screen.
        return "\n".join([self.header_line, self.sequence_line, self.separator_line, self.quality_line])
        # returns the header line, sequence line, separator line and quality line of the read on separate lines


class ReadPair:
    """Sequenced read mates and associated data."""

    def __init__(self, forward_read, reverse_read):
        """Constructor for ReadPair"""
        self.forward_read = forward_read
        self.reverse_read = reverse_read
        # If the read was assigned to a sample then set to the sample_id, otherwise set to False
        # Default to False (unassigned read)
        self.sample = False

    def trim(self, start, stop):
        """Crops the read sequence and associated quality_line keeping positions start to stop.
        Be careful to give the positions in Python position numbering (starts at 0). The pain of calculating the proper
        coordinates before calling the function will save a lot of time avoiding re-calculating them every time the
        function is called.

        Args:
            self, start, stop

        Returns:
            None
        """
        self.forward_read.trim(start, stop)
        self.reverse_read.trim(start, stop)

    def define_quality_status_forward(self, threshold , percentage):
        """Sets the quality_status attribute of the forward read according to whether stricly more than percentage of
        bases in the read are strictly below phred.
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
        self.forward_read.define_quality_status_forward(threshold, percentage)

    def define_quality_status_reverse(self, threshold , percentage):
        """Sets the quality_status attribute of the forward read according to whether stricly more than percentage of
        bases in the read are strictly below phred.
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
        self.reverse_read.define_quality_status_forward(threshold, percentage)

    def define_adapter_presence_substitutions_only_forward(self, adapter, max_substitutions):
        """Sets the adapter_absent attribute according to whether a match is found with a number of substitutions
        less or equal to
        edit_threshold.

        :rtype : None
        Args:
            self, adapter, edit_threshold

        Returns:
            None

        Args:
            self, adapter, max_substitutions

        Returns:
            None
        """
        self.forward_read.define_adapter_presence_substitutions_only(adapter, max_substitutions)

    def define_adapter_presence_substitutions_only_reverse(self, adapter, max_substitutions):
        """Sets the adapter_absent attribute according to whether a match is found with a number of substitutions
        less or equal to
        edit_threshold.

        :rtype : None
        Args:
            self, adapter, edit_threshold

        Returns:
            None

        Args:
            self, adapter, max_substitutions

        Returns:
            None
        """
        self.reverse_read.define_adapter_presence_substitutions_only(adapter, max_substitutions)

    def __str__(self):
        """String representation of the read pair.

        Args:


        Returns:
            None
        """
        return "%s\n%s"% (self.forward_read, self.reverse_read)


class Adapter:
    """Adapter sequence and associated data"""

    def __init__(self, identifier, sequence, mate):
        """Constructor for Adapter"""
        self.identifier = identifier
        self.sequence = sequence
        self.mate = mate

    def __str__(self):
        """This function defines a string representation displayed when calling the print function on an Adapter
        object"""
        # Uses the tab separator to separate the different elements of the Adapter
        return '\t'.join([self.sequence, self.identifier, self.mate])
        # returns the header line, sequence line, separator line and quality line of the read on separate lines

