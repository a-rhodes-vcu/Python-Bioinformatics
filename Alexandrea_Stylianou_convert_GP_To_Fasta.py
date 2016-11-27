from __future__ import print_function

import re
import textwrap

class Question1(object):
    a = re.compile("DEFINITION\s(\W.+)(\[|,)")
    b = re.compile("SOURCE\s+(\w.*)")
    c = re.compile("VERSION\s+(\w+\W\w)\s+GI\W(\w+)")
    d = re.compile("^\d{1,3}\s(\w.*)")

    def __init__(self, file):
        self.file = file
        self.header = []
        self.seq_string = '' #i was only able print my sequence string out with self.seq_string.
                             #this took FOREVER TO FIGURE OUT

    def parse_header(self):
        first_def1 = source_line = gi_1 = gi_2 = seq_string = None #setting all my variables to None
        genplet_data = None #this is my tuple, set that to None
        for lines in self.file:
            lines = lines.strip() #strip out white space
            if lines.startswith("LOCUS"):
                if first_def1 and source_line and gi_1 and gi_2 and seq_string: #took this from blast parser.
                    genplet_data = (gi_2, gi_1, first_def1, source_line,seq_string)
                    self.header.append(Alignments(*genplet_data))
                    self.seq_string = '' # this is a flag, once the line is empty stop
            elif lines.startswith("DEFINITION"): #if line starts with definition
                accession_id = self.a.search(lines)
                if accession_id: #if the regular expression exisits
                    first_def1 = accession_id.group(1) #grab group 1, and so forth
                    #remove_comma = re.search("(\w.*)(,.*)",first_def) #trying to remove pesky comma.
                    #if remove_comma:                                  #gave up
                        #first_def1 = remove_comma.group(1)
            elif lines.startswith("SOURCE"):
                source_info = self.b.search(lines)
                if source_info:
                    source_line = source_info.group(1)
            elif lines.startswith("VERSION"):
                version_info = self.c.search(lines)
                if version_info:
                    gi_1 = version_info.group(1)
                    gi_2 = version_info.group(2)
            elif re.search("^\d+",lines): #if line starts with a digit
                seq_match = self.d.search(lines)
                if seq_match:
                    #first attempt was this, which was not working. ended up putting everything into one big line.
                    #uppercase_seq = full_seq.upper()
                    #remove_space = uppercase_seq.replace(" ", "")
                    #seq_string += remove_space
                    self.seq_string += seq_match.group(1) #add each sequence to self.seq_string
                    seq_string = self.seq_string.replace(" ","").upper() #remove the space, change to uppercase
        genplet_data = (gi_2, gi_1, first_def1, source_line, seq_string)
        self.header.append(Alignments(*genplet_data))
        return self.header

    def write_report(self, header, outfile): #call this method when writing report to outfile
        with open(outfile, 'w') as ofh:
            for a in header:
                print(a, file=ofh, end='')


class Alignments(object):
    def __init__(self,gi_2,gi_1,first_def1,source_line,seq_string):
        self.gi_2 = gi_2
        self._gi_1 = gi_1
        self.first_def1 = first_def1
        self.source_line = source_line
        self.seq_string = seq_string

    def __str__(self):
        return ">gi|{}|gb|{}|{}\n[{}]\n{}\n".format(self.gi_2,
                                                    self._gi_1,self.first_def1,self.source_line,
                                                    textwrap.fill(self.seq_string, 70)) #70 residues per line

def main():
    file_name = "O104_H4_GP.txt" #please add file name here
    with open(file_name) as infile: #took all of this from blast parse
        report = Question1(infile)
        alignments = report.parse_header()
        report.write_report(alignments, "0104_H4_FASTA.txt")
        #for thing in alignments: #only need this if printing to console
            #print(thing)

if __name__ == "__main__":
    main()

