__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module os allows to access OS-specific line separator.
import os

class Read:
    """Sequenced read and associated data"""

    def __init__(self, header_line, sequence_line, separator_line, quality_line):
        """Constructor for Read"""
        self.header_line = header_line
        self.sequence_line = sequence_line
        self.separator_line = separator_line
        self.quality_line = quality_line

    def trim(self, five_start, three_end):
        """Crops the read sequence and associated quality_line keeping positions five_start to three_end.
        Be careful to give the positions in Python position numbering (starts at 0).
        """
        self.sequence_line = self.sequence_line[five_start:three_end]
        self.quality_line = self.quality_line[five_start:three_end]

    def __str__(self):
        """This function defines the string representation displayed when calling the print function on a Read object"""
        return os.linesep.join([self.header_line, self.sequence_line, self.separator_line, self.quality_line])

