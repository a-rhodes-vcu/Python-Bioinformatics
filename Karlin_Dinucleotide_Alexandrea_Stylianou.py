from __future__ import division
from FastA_V2 import FastA
import matplotlib.pyplot as plt


class Karlin_dinucleotide(object):

    """The complete DNA sequence for Vibrio cholerae O1 biovar El Tor str. N16961 chromosome II can be found at:
    https://www.ncbi.nlm.nih.gov/nuccore/15640032?report=fasta

    This program should calculate local and genome-wide nucleotide and dinucleotide frequencies, then using the
    method of Karlin, calculate average absolute dinucleotide relative abundance differences between relative
    dinucleotide abundances found in each window vs. the relative dinucleotide abundances for the genome as a whole.
    I've given you lots of clues as to the methods that are required and often even to the variables names you
    might wish to use. Let's play fill in the blanks!
    """

    def __init__(self, sequence_file='Vib_Cholerae_Fasta.fasta', win_size=20000, skip=1000):
        """
        The initializer accepts three arguments. The first is the file name, the second is the desired
        """

        self.ORDER = 2  # I am capitalizing this to make it clear it is a constant.

        deltas = []  # A list where will store the delta values associated with each window.


        # A tip off to some of the variables you might wish to use later...
        genome_NT_freqs = {}
        genome_DN_freqs = {}
        window_NT_freqs = {}
        window_DN_freqs = {}

        # First we will appropriately count all the dinucleotides in the whole genome, assuming double-stranded DNA

        fast_A_object = FastA(sequence_file)

        annotation, sequence = fast_A_object.next()  # This time we're expecting a fast A file with only one sequence!
        # instead of iterating when we know we only want one thing, we'll just ask for it directly with .next()

        # The Vibrio sequence contains a few characters from the nucleotide ambiguity code. Cheat and get rid of these!
        # See http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
        sequence = sequence.replace('Y', 'C')
        sequence = sequence.replace('R', 'A')
        sequence = sequence.replace('K', 'T')
        sequence = sequence.replace('N', 'T')
        sequence = sequence.replace('W', 'A')
        sequence = sequence.replace('S', 'G')
        sequence = sequence.replace('M', 'C')

        # Calculate GC frequencies and dinucleotide frequencies for the whole genome!
        genome_NT_freqs, genome_DN_freqs = self._Calc_NT_and_DN_frequencies(sequence)

        print "Genome-wide nucleotide frequencies are:\n", genome_NT_freqs
        print "Sum of genome-wide nucleotide frequencies is", sum(genome_NT_freqs.values())
        print "Genome-wide dinculeotide frequencies are:\n", genome_DN_freqs
        print "Sum of genome-wide dinucleotide frequencies is", sum(genome_DN_freqs.values())

        # Once we have these, we can calculate the dinucleotide relative abundance
        genome_DN_abundance = self._Calculate_DN_relative_abundance(genome_NT_freqs, genome_DN_freqs)

        print "Genome-wide relative dinucleotide abundances are:\n", genome_DN_abundance
        print "Average of genome-wide relative dinucleotide abundances (should be very close to 1.00) are:\n", \
            sum(genome_DN_abundance.values()) / 16

        # Now that we have what we need for the whole genome, we can start busting it into chunks and calculating local
        # GC percentage and dinucleotide frequencies, as well as the relative dinucleotide abundances
        # for each of many windows.  These can further be used to calculate the differences
        GC_freq = []

        for i in xrange(0, len(sequence) - win_size + 1, skip):
            print "Working on window", int(i / skip)
            # First we get the GC and DN frequencies for the window

            GC = 0
            c_value = 0
            g_value = 0
            window_sequence = sequence[i : i + int(win_size)]
            window_NT_freqs,window_DN_freqs = self._Calc_NT_and_DN_frequencies(window_sequence)
            ############################### grab the frequencies collected from the window dict    ######################
            GC=float(window_NT_freqs["G"] + window_NT_freqs["C"]) *100 # do the calculations

            GC_freq.append(GC)
            print("This is the GC frequency",GC_freq)
            ############################# put GC frequencies into a list    ##################################


            window_DN_abundance = self._Calculate_DN_relative_abundance(window_NT_freqs, window_DN_freqs)

            # Finally, we can calculate the local delta between the genome and window according to the formula from
            # Karlin.  Remember, this is really just an average of the absolute differences between the relative
            # dinucleotide abundances of the window vs. relative dinucleotide abundances of the genome as a whole.
            window_delta = 0
            for k,v in window_DN_abundance.items():
                window_delta += abs(v-genome_DN_abundance[k]) #using Karlin's formula, the relative abundance distance is made.
            window_delta = (1/16) * window_delta

            #Since we are going to want to plot these values later, you should probably stash each of the values that
            # you calculate into a list of deltas that we can use later with matplotlib

            deltas.append(window_delta)
            #print "this is the deltas",deltas

        # Once we are finished with the loop, we can print out or plot the deltas.  Here is some code that uses
        # matplotlib to make the charts like in Karlin's paper... If you have a list with all the deltas, this
        # should just work.
        ############################# di nuc graph
        x_values = [i * skip for i in xrange(len(deltas))]
        plt.plot(x_values, deltas)
        plt.ylabel('local delta relative to genome delta')
        plt.xlabel('Genome position')
        plt.show()
        ##########################  GC graph
        x_values = [i * skip for i in xrange(len(GC_freq))]
        plt.plot(x_values, GC_freq)
        plt.ylabel('GC % content')
        plt.xlabel('Genome position')
        plt.show()

    def _Calc_NT_and_DN_frequencies(self, sequence):
        """
        This method takes as an argument a sequence of arbitrary length that is presumed to be double stranded,
        and should return a tuple of two dicts: the first dict should have keys that are nucleotides and values that are
        nucleotide frequencies, while second should have keys that dinucleotides and values that are dinucleotides.
        """

        # Here are some dicts that might prove handy
        dinuc_RC = {'CG': 'CG', 'GC': 'GC', 'TA': 'TA', 'AT': 'AT', 'CC': 'GG', 'GG': 'CC', 'TT': 'AA', 'AA': 'TT',
                    'TG': 'CA', 'CA': 'TG', 'AG': 'CT', 'CT': 'AG', 'AC': 'GT', 'GT': 'AC', 'GA': 'TC', 'TC': 'GA'}

        nuc_RC = {'A': 'T', 'C': 'G', 'G' : 'C', 'T': 'A'}

        # These are the dicts that will be populated with raw counts and frequencies for dinucleotides and NTs
        DN_counts = {}
        NT_counts = {}

        DN_freq_dict = {}
        NT_freq_dict = {}

        di_nuc_denom = (len(sequence) - self.ORDER + 1) *2 #this is the effective length, went over this in class

        for i in xrange(len(sequence) - self.ORDER + 1):
            #getting 5-3 strand
            direct_word = sequence[i: i + self.ORDER]
            #getting 3-5 strand
            revcom_word = dinuc_RC[direct_word]

            #get single nucleotides
            direct_NT = direct_word[0]
            # reverse single nucleotides
            revcom_NT = nuc_RC[direct_NT]

            #getting 5-3 strand
            direct_DC = direct_word
            #getting 3-5 strand
            revcom_DC = dinuc_RC[direct_DC]

            try:
                #do the counts
                NT_counts[direct_NT] += 1
                NT_counts[revcom_NT] += 1
                #################### get the dinuc counts
                DN_counts[direct_DC] += 1
                DN_counts[revcom_DC] += 1

            except KeyError:
                #if doesn't exisit, then create a key, set to one
                NT_counts[direct_NT] = 1
                NT_counts[revcom_NT] = 1

                DN_counts[direct_DC] = 1
                DN_counts[revcom_DC] = 1

        # To get the count just right, once the loop is over you may need to add the last nucleotide to the count..
        NT_counts[direct_word[1]] += 1
        NT_counts[revcom_word[1]] += 1

      #now we get the single nucleotide freq. we have the counts stored in a dic. then divide by the length times 2
        for k,v in NT_counts.items():
            NT_freq_dict[k] = float(v/(len(sequence)*2)) #for this multiply the length by two, to take into account the 3-5 strand
        print("nt freq dict",NT_freq_dict)

        #create a new dictionary with di nuc frequencies. we  have the counts, now divide by the effective length
        for k,v in DN_counts.items():
            DN_freq_dict[k] = float(v/(di_nuc_denom)) #divide the values by the effective length

        return NT_freq_dict, DN_freq_dict

    def _Calculate_DN_relative_abundance(self, NT_freq_dict, DN_freq_dict):
        """
        This method accepts as arguments two dicts, the first which describes single nucleotide frequencies,
        and the second describing dincucleotide frequencies.  The method should then calculate a relative dinucleotdide
        abundance using the formula described by Karlin.  Remember, this is essentially the ratio of the Observed
        over the Expected number of dinucleotides, where the Expected number is how many we should expect to have
        observed under our neutral model that individual nucletides were randomly (i.e. independently) assorting
        with one another.. we can also think of this expected value as a zeroth order Markov model.  The keys and values
        are as described in the Calc_GC_and_DN_frequencies method.  The return value here should of course be a dict
        associating each dinucleotide (key) with the relative dinucleotide abundance (value) for that dinucleotide.
        In principle a clever programmer could reduce this whole problem to one line of code with a return statement
        and a dict comprehension.  Are you up for it?

        """
        ################### observed over expected  #######################
        DN_Dict = {k: float(v) / (NT_freq_dict[k[0]] * NT_freq_dict[k[1]]) for k,v in DN_freq_dict.items()}
        return DN_Dict

A = Karlin_dinucleotide('Vib_Cholerae_Fasta.fasta', win_size=50000, skip=5000)