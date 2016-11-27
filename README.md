# Python-Bioinformatics

Hello World!

This document is designed to give a brief description of the programs uploaded, and what files they use. These programs were either written by me, or modified by me.

1.Alexandrea_Stylianou_convert_GP_to_Fasta.py:
 
 converts a gen bank file: https://www.ncbi.nlm.nih.gov/protein/559119974?report=genpept to a fasta file (shown below):
>AA 1 
MDSLLFDRRPDIQLPPQSC*HKYPRKSERETD*LNTGSDRKVSV*HIYV*RSHG*VIAKLMQQR*NQRILLGFYYCNGSMYQEP*SIGRQVFYLIFVSRRATRIYLHVR*RIQAT*MFYRLSPPSSKFFKCVARRVHHRLNSVS*MPHGAWMPLYEFSTVIQKTVYLLNKCAGLSLLVTPLLTR**GDLKCRNWNIHFRQS**CSKQTTDQLRYVY*PPEFLMWSIPTDSFTLSGLIKCDTYGFPLQTFPNINKQ*FAFSMCSTKR*AQIPIQEDNKGDSYHRGRHS*FYYRGEE*VTDIHMHGMILTTSY*WPHVQGCSSVGRTLCTP*IVV

2.Alexandrea_Stylianou_bias_analyzer: 
Creates 6000 simulated coding regions, 999 nucleotides long, then compare the codon frequency to the codon frequency found in     cyanobacteria. 
  
3.Stylianou_BLAST_Prot:
This program was modified to simulate Protein Blast. In essence, this program utilizes a PAM120 table, and scores protein comparisons     based on an exact, good, poor, or bad match. 

4.Stylianou_Linkage_Cluster:
Linkage clustering is a modified program that uses Pearson Coefficent in order perform single linkage clustering or complete linkage       clustering. This program uses a text file of microarray data.  
