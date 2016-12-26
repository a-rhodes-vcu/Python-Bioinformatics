from pickle import load
import AA_module
import FastA




class ThreadProtein(object):

    """ThreadProtein
    VERSION:    Version 3.0 (3 Nov 2015) by Paul Fawcett, but very closely modelled after an earlier PERL program by
                Jeff Elhai

    PURPOSE:    Threads the amino acid sequence of one protein
                through the three-dimensional structure of another

    INPUT FILES:

        Structure file: PDB format. Must have only one chain, corresponding to known sequence.

        Alignment file: FastA format (output as PIR format by Clustal)

            Leading and lagging blanks are OK

            SEQUENCES MUST BE ALIGNED! (i.e. with gaps if necessary)
            First sequence must be that also represented by structure file = Known sequence
            Second sequence must be that to be superimposed on known structure = threaded sequence

            > [header line for sequence with known structure]
            [optional blank line]
            [amino acid sequence, with gaps (-) for alignment]
            * [optional]
            > [header line for sequence to be threaded]
            [optional blank line]
            [amino acid sequence, with gaps (-) for alignment]
            * [optional]

    OUTPUT FILES: Description of format of output files

        New structure file: PDB format

            NO CHANGE: Amino acids common between the two sequences go in chain A
            REPLACEMENTS: Amino acids changed in threaded sequence go in chain B
                Amino acid called new name in first atom of residue (others unchanged)
            INSERTIONS: Amino acids in threaded sequence but not original go in chain C
                No attempt made to find positions of the new atoms.
                Only first atom given
                Inserted atom given meaningless number
            DELETIONS: Amino acids in original but not in threaded sequence go in chain D

            Informational lines (e.g. Title, source, etc) copied unchanged
            therefore information may no longer be appropriate!!
    """

    def __init__(self, structure_file='1ei1.pdb', alignment_file='Consensus_Ecoli_Alignment.pir', output_file='Consensus_Threading.pdb',list_conserved_numbers="info.txt",threshold_number=3.0,highly_conserved_threshold = 4.0):

        self.structure_file = structure_file
        self.alignment_file = alignment_file
        self.output_file = output_file
        self.threshold_number= threshold_number
        self.highly_conserved_threshold = highly_conserved_threshold
        self.translation_dict = {}


        headers, sequences = FastA.Read_FastA_sequences(alignment_file)

        self.known_seq_header = headers[0]
        self.threaded_seq_header = headers[0]


        self.known_seq = sequences[0]
        self.threaded_seq = sequences[1]

        self.information = []
        with open(list_conserved_numbers, 'r') as f:
            self.information = load(f)
        self.comparison_summary = self.__Analyze_alignment()
        comparison_summary_list = []
        for i in self.comparison_summary:
            comparison_summary_list.append(i)
        self.comparison_summary2 = self.__Modify_comparison_summary(comparison_summary_list)
        print(self.comparison_summary2)

        self.__Create_new_structure_file()
        print "Finished writing new PDB file", self.output_file



    def __Analyze_alignment(self):

        """ Analyze alignment at each position
        Put I (insertion) where gap in sequence with known structure
        Put D (delection) where gap in sequence to be threaded
        Put M (match) where both sequences match
        Put R (replacement) where sequences differ

        relies on instance variable self.known_seq and self.threaded_seq
        """

        comparison_summary = ''
        gap_char = '-'

        print "\nAnalyzing alignment:\n"
        print self.known_seq
        print self.threaded_seq


        if not len(self.known_seq) == len(self.threaded_seq):

            quit("The two sequences must be the same length. Don't forget to align them!")

        print "Length of input sequences is", len(self.known_seq)

        ## First we create a string that summarizes the alignment in a position-by-position fashion

        for i in xrange(len(self.known_seq)):

            known_char = self.known_seq[i]
            threaded_char = self.threaded_seq[i]

            if known_char == gap_char:
                comparison_summary += "I"
            elif threaded_char == gap_char:
                comparison_summary += "D"
            elif known_char == threaded_char:
                comparison_summary += "M"
            else:
                comparison_summary += "R"

        comparison_summary += "M"           # Adds a final match so as to correctly handle the termination sequence


        return comparison_summary

    def __Modify_comparison_summary(self,comp_list):
        conv_position = 0

        #did this in class, populate dictioanry: 0: 53, 1: 54, 2: 55, 3: 56, 4: 57, 5: 58...

        for summary_position, letter in enumerate(self.comparison_summary):
            if letter != 'D': #collect the counts of all non-d's
                conv_position += 1
            self.translation_dict[conv_position] = summary_position #make a dictionary with corresponding positions
        #print (self.translation_dict)

        for i,v in enumerate(self.information):
            if v>self.threshold_number: #if v is over a threshold then
                comp_list[self.translation_dict[i]+1] = '^' #have to move over indexes by one in order to get star in correct place

        for i,v in enumerate(self.information):
            if v>self.highly_conserved_threshold: #if v is over a threshold then
                comp_list[self.translation_dict[i]+1] = '+'

        mutation_comparison_summary = ''.join(comp_list) #turn list back into a string
        comparison_summary = mutation_comparison_summary

        return comparison_summary

    def __Create_new_structure_file(self):


        gap_char = '-'
        gaps_reached_in_known_seq = 0
        terminus_found = False
        chains_printed = False
        previous_line = ''

        last_residue = len(self.known_seq.translate(None, gap_char))
        atom_lines = dict()

        atom_lines['A'] = []
        atom_lines['B'] = []
        atom_lines['C'] = []
        atom_lines['D'] = []
        atom_lines['E'] = [] #fifth chain, to hold mutant information in
        atom_lines['F'] = []

        with open(self.structure_file, 'rU') as old_structure: #opening the structur file
            with open(self.output_file, 'w') as new_structure: #opening new output file. spit out a modified version

                line = old_structure.readline()
                while line:

                    line = line.strip()
                    line_type = line[0:4]

                    if line_type == 'ATOM' or line_type == 'TER ': #if line time is atom, or ter then pick up a reside from positin on residue from 22 to 26. relying on structure of PDB file. hard coded.
                        residue = int(line[22:26]) - 1 + gaps_reached_in_known_seq

                        if line_type == 'TER ':     # PDB file numbers terminal residue
                            residue += 1            # same as last residue. Fix this. if a match nothing is changed, rights back out.
                            terminus_found = True   #match, rewriting old line, nothing has changed

                        #print residue, "is", self.comparison_summary[residue]
                        residue_status = self.comparison_summary2[residue]

                        if residue_status == 'M':

                            print >>new_structure, line
                            previous_line = line
                            line = old_structure.readline()

                        elif residue_status == 'R':

                            new_residue = self.threaded_seq[residue]
                            #print(self.threaded_seq[residue])
                            #print(line)
                            line = self.__Update_residue_name(line, new_residue)
                            line = self.__Update_chain_ID(line, 'B')
                            atom_lines['B'].append(line)
                            previous_line = line
                            line = old_structure.readline()

                        elif residue_status == '^': #this is the new chain, for bit values above 3
                            newline = self.__Update_chain_ID(line, 'E')
                            atom_lines['E'].append(newline)
                            previous_line = line
                            line = old_structure.readline()

                        elif residue_status == '+': #this is the new chain, for bit values above 4
                            newline = self.__Update_chain_ID(line, 'F')
                            atom_lines['F'].append(newline)
                            previous_line = line
                            line = old_structure.readline()


                        elif residue_status == 'D': #we had a deletion with respect to strep. yellow amino acids don't exisit  in m loti. strep is our structure.
                            newline = self.__Update_chain_ID(line, 'D')
                            atom_lines['D'].append(newline)
                            previous_line = line
                            line = old_structure.readline()

                        elif residue_status == 'I': #not sure what to do, exisit in loti, but not in strep.
                            last_residue += 1
                            new_residue = self.threaded_seq[residue]
                            if line_type == 'TER ':
                                newline = previous_line
                            else:
                                newline = line

                            newline = self.__Update_residue_name(newline, new_residue)
                            newline = self.__Update_chain_ID(newline, 'C')
                            newline = self.__Update_residue_number(newline, last_residue)
                            gaps_reached_in_known_seq += 1
                            atom_lines['C'].append(newline)

                        else:
                            print 'Error in:\n', line, '\n\tResidue', residue, '\n\tResidue status:', residue_status

                            line = old_structure.readline()

                    else:
                        if terminus_found and not chains_printed:

                            self.__Output_extra_chains(atom_lines, new_structure)
                            chains_printed = True
                        print >>new_structure, line
                        line = old_structure.readline()

    def __Update_residue_name(self, line, residue_name):

        line = list(line)       # string type doesn't support item assignment, so convert to list

        if residue_name == 'INS' or residue_name == 'DEL':
            line[17:20] = residue_name
        else:
            line[17:20] = AA_module.One_letter_to_three_letter_code(residue_name)

        return ''.join(line)
        #return str(line) #something is being cast to a string.

    def __Update_chain_ID(self, line, chainID):

        line = list(line)       # string type doesn't support item assignment, so convert to list

        if not len(chainID) == 1:  # Chain ID should ALWAYS be length 1
            print "ERROR: Chain ID", chainID, "not of length 1"
            quit()
        line[21] = chainID #overright one positon in the list, to something new

        return ''.join(line)

    def __Update_residue_number(self, line, residue_number):

        line = list(line)       # string type doesn't support item assignment, so convert to list

        new_number = "%4d" % residue_number
        line[22:26] = new_number
        return ''.join(line)

    def __Output_extra_chains(self, atom_lines, filehandle):

        #print "Outputing remaining chains"

        for chain in sorted(atom_lines.keys()):

            for atom in xrange(len(atom_lines[chain]) - 1):

                print >>filehandle, atom_lines[chain][atom]


def main():

    ThreadProtein()  # Just use the default filenames for now

if __name__ == '__main__':

    main()