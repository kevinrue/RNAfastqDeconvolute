__author__ = 'David Magee, Carolina Correia, and Kevin Rue-Albrecht'
__copyright__ = "Copyright 2014, GPLv2"

"""Empty docstring
"""

# Module gzip allows to parse compressed files without having to first decompress them.
import gzip
# Custom module to store the data of a read
from SeqDataTypes import Read


class FastqgzParser:
    """Parses compressed fastq files."""

    def __init__(self, filename):
        """Constructor for FastqqzParser
        :rtype : A parser of class FastqgzParser
        """
        self.filename = filename
        self.file_handler = gzip.open(self.filename, 'rb')

    def open(self):
        """Opens the file stream.

        Args:
            self

        Returns:
            None
        """
        self.file_handler = gzip.open(self.filename, "r")

    def nextRead(self):
        """Parses and return the next sequenced read in the file.

        Args:
            self

        Returns:
            The next sequenced read in the file.
        """
        header_line = self.file_handler.readline().strip()
        sequence_line = self.file_handler.readline().strip()
        separator_lne = self.file_handler.readline().strip()
        quality_line = self.file_handler.readline().strip()
        return Read.Read(header_line, sequence_line, separator_lne, quality_line)

    def close(self):
        """Closes the file stream.
        
        Args:
            self
           
        Returns:
            None
        """
        self.file_handler.close()