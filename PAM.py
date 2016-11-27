from __future__ import division
from math import log10


class PAM(object):

    """Following in the footsteps of Margaret Dayhoff, we create a PAM matrices of arbitrary evolutionary
    distance and create a dict useful for scoring alignments in the Smith-Waterman or Blast algorithms
    """

    def __init__(self, N=250, filename='PAM1.txt', aa_frequencies=None):

        self.PAM = {}
        # The scoring dict we are going to create and print out. Keys are in the form (aa1, aa2) in the one-letter code
        # values can be directly used for scoring in either a protein Smith-Waterman or PBLAST type algorithm

        self.N = N                    # The evolutionary distance of the PAM matrix we desire
        self.filename = filename      # We skip quite a few steps, and start with a PAM1 probability matrix

        if aa_frequencies is None:

            self.aa_frequencies = {"G": 0.089, "R": 0.041, "A": 0.087, "N": 0.040, "L": 0.085, "F": 0.040, "K": 0.081,
                                   "Q": 0.038, "S": 0.070, "I": 0.037, "V": 0.065, "H": 0.034, "T": 0.058, "C": 0.033,
                                   "P": 0.051, "Y": 0.030, "E": 0.050, "M": 0.015, "D": 0.047, "W": 0.010}
        else:

            self.aa_frequencies = aa_frequencies  # OK, so the normalized amino acid frequencies are a gimme too!

    def Build_PAMN(self):

        PAM1, alphabet = self.__read_PAM1()       # Read in the PAM1 file

        size = len(alphabet)  # Record how many aa were in the input file
        # Some PAM table include ambiguity characters, so we cannot assume 20

        PAM1 = [[element / 10000 for element in row] for row in PAM1]
        # This sneaky-looking nested list comprehension just divides each element in our 2D matrix by 10000
        # nested loop comprehensions are a compact and idiomatically Pythonic way to perform some simple operation on
        # each element of a two (or more) dimensional list

        PAMN = list(PAM1)       # We need to start off by multiplying PAM1 by itself, so we need a copy

        # Note that the above gives a true "clone" of PAM1, not just the "shallow copy" we would using just PAMN = PAM1
        # Do some experimenting and discover the issues with making a shallow copy of a list.  For instance, try
        # appending a new item to the shallow copy. Print out the original to see the (perhaps unanticipated) effect

        for i in range(self.N - 1):  # Repeatedly multiply by PAM1 to end up with a true PAMN table
                                      # We started already with PAM1, so only need to multiply N - 1 times

            PAMN = self.__MatrixMultiply(PAMN, PAM1)

        for i in range(size):  # Convert to a log odds formulation. Could you redo this as a list comprehension?
            for j in range(size):
                PAMN[i][j] = log10(PAMN[i][j] / self.aa_frequencies[alphabet[i]])

        for i in range(size):  # Average the reciprocal values, multiply by 10, round, and build the final dict
            for j in range(size):
                self.PAM[alphabet[i], alphabet[j]] = int(round(10 * ((PAMN[i][j] + PAMN[j][i]) / 2)))

        return self.PAM

    def __read_PAM1(self):

        """Read in the raw PAM1 file into a two-dimensional list"""

        PAM1 = []       # List for the starting PAM1 matrix as read from a file

        with open(self.filename, 'rU') as pam1_file:

            line = pam1_file.readline()  # read the header line of the PAM file indicating the aa ordering in the table
            line = line.strip()          # chew off the endline character

            alphabet = line.split('\t')
            # Break up the tab-delimited line into an array indicating the order the amino acids appear in the table

            del alphabet[0]   # Throws away the first element, which is not an amino acid, just a label.
            # alphabet will reflect the order of the aa as they appeared in the input file

            line = pam1_file.readline()    # Get the next line, which is first line of matrix

            while line:

                line = line.strip()                             # Strip off the trailing line feed
                line = line.split('\t')                         # Grab the whole tab-delimited line as an array
                del line[0]                                     # Dump the first element, which is an aa, not a number
                line = [int(element) for element in line]       # Convert the entries from strings to numbers
                PAM1.append(line)
                line = pam1_file.readline()

        return PAM1, alphabet

    def PAM_dump(self):

        """Print out the final scoring dict we have created. """

        ordered_symbols = self.aa_frequencies.keys()        # If I was a bit less lazy I would get this from alphabet
        #ordered_symbols.sorted()
        ordered_symbols = sorted(ordered_symbols)

        for aa in ordered_symbols:
            print ("\t", {}).format(aa)
        print ('\n')
        for aa1 in ordered_symbols:
            print ({},'\t'.format(aa1))
            for aa2 in ordered_symbols:
                print ({}, '\t').format(self.PAM[aa1, aa2])
            print ('\n')

    def __MatrixMultiply(self, arrayA, arrayB):

        """Multiply two matrices, A and B, and return matrix C. May dot products haunt your dreams!"""

        rowsA = len(arrayA)         # Rows of A
        colsB = len(arrayB[0])      # Columns of B

        if rowsA != colsB:

            raise ValueError("math domain error")

        new_array = [[0 for i in range(rowsA)] for j in range(colsB)]  # initialize an empty array of correct size

        rowcols = len(arrayA[0])                         # Number of columns in A and rows in B

        for i in range(rowsA):                          # iterate over rows of A
                for j in range(colsB):                  # iterate over columns of B.
                        for k in range(rowcols):        # iterate over columns of A and rows of B

                                new_array[i][j] += arrayA[i][k] * arrayB[k][j]

        return new_array


def main():

    A = PAM(N=200)
    A.Build_PAMN()
    A.PAM_dump()

if __name__ == '__main__':

    main()