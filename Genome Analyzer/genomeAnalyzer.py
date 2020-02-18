#!/usr/bin/env python3
# Name: Justin Juang (jjuang)
# Group Members: None

'''
genomeAnalyzer.py program imports classes NucParams and FastAreader from sequencAnalysis module, 
calculates and prints out the sequence length, G and C nucleotide content, and relative codon 
frequency and codon count for each codon from a fastA file.

Sequence length: total length of genome sequence (sum of all nucleotides). 
GC content: number of G + C nucleotides divided by the total number of all nucleotides.
Codon fequency: number of particular codon divided the number of amino acid of which the codon codes for.
'''

from sequenceAnalysis import NucParams, FastAreader

def main ():
    '''
    Function main that uses imported classes NucParams and FastAreader to print out genome length, GC
    nucleotide content, and codon frequency for each codon. 
    '''

    myReader = FastAreader() #instance of fastA file 
    myNuc = NucParams() #instance of NucParams
    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)

    #Prints sequence length.
    seqLength = myNuc.nucCount()
    mbSeqLength = seqLength/1000000
    print ("sequence length: {0:0.2f}Mb".format(mbSeqLength),"\n")

    #Printes GC content.
    nucComp = myNuc.nucComposition()
    gAndC = nucComp["G"] + nucComp["C"]
    gcContent = gAndC/myNuc.nucCount()
    print ("GC content = {0:0.1%}%".format(gcContent),"\n")

    #Prints relative codon frequency and codon count for each codon.
    aaComp = myNuc.aaComposition() 
    rnaCodonComp = myNuc.codonComposition()
    sortedRNATable = sorted(myNuc.rnaCodonTable.items(), key=lambda x:x[1]) #sorts aa by alphabetical order.  

    for codon, aa in sortedRNATable: #for each codon in rna codon table
        aaNumber = aaComp[myNuc.rnaCodonTable[codon]] #denominator for relative codon frequency calculation.
        codonNumber = myNuc.rnaCodonCompDict[codon] #numerator 
        if aaNumber > 0: 
            codonFreq = codonNumber/aaNumber #relative codon frequency
        elif aaNumber == 0:
            codonFreq = 0
        print ("{:s} : {:s} {:5.1f} ({:6d})".format(codon, aa, codonFreq * 100, codonNumber), end = " ")

if __name__ == "__main__":
    main()
    
