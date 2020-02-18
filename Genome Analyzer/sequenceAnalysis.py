#!/usr/bin/env python3
# Name: jjuang (jjuang)
# Group Members: None

'''
The program reads and calculates the physical-chemical properties of a user input
protein sequence. The program will read in the protein sequence and print out the
following: number of amino acids and total molecular weight, molar extinction 
coefficient and mass extinction coefficient, theoretical isoelectric point (pI), 
and amino acid composition.

Example:
    Input: VLSPADKTNVKAAW
    Outpur: Number of Amino Acids: 14
            Molecular Weight: 1499.7
            molar Extinction coefficient: 5500.00
            mass Extinction coefficient: 3.67
            Theoretical pI: 9.88
            Amino acid composition (in seperate lines): A = 21.43% C = 0.00% D = 7.14% 
            E = 0.00% F = 0.00% G = 0.00% H = 0.00% I = 0.00% K = 14.29% L = 7.14% M = 0.00% 
            N = 7.14% P = 7.14% Q = 0.00% R = 0.00% S = 7.14% T = 7.14%  V = 14.29% W = 7.14% 
            Y = 0.00%
'''

class NucParams:
    rnaCodonTable = {
                    # RNA codon table
                    # U
                    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
                    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
                    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
                    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
                    # C
                    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
                    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
                    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
                    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
                    # A
                    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
                    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
                    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
                    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
                    # G
                    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
                    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
                    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
                    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
                    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()} #dna codon table

    def __init__ (self, inString=''): #aa dict, nucleotides dic, and codon dict
        self.aaCompDict = {
                            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0, 
                            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
                            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
                            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0
                            }
       
        self.nucCompDict = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'U': 0, 'N':0}

        self.rnaCodonCompDict = {
                    'UUU': 0, 'UCU': 0, 'UAU': 0, 'UGU': 0,  
                    'UUC': 0, 'UCC': 0, 'UAC': 0, 'UGC': 0,  
                    'UUA': 0, 'UCA': 0, 'UAA': 0, 'UGA': 0,  
                    'UUG': 0, 'UCG': 0, 'UAG': 0, 'UGG': 0,  
                    
                    'CUU': 0, 'CCU': 0, 'CAU': 0, 'CGU': 0,  
                    'CUC': 0, 'CCC': 0, 'CAC': 0, 'CGC': 0,  
                    'CUA': 0, 'CCA': 0, 'CAA': 0, 'CGA': 0,  
                    'CUG': 0, 'CCG': 0, 'CAG': 0, 'CGG': 0,  
                    
                    'AUU': 0, 'ACU': 0, 'AAU': 0, 'AGU': 0, 
                    'AUC': 0, 'ACC': 0, 'AAC': 0, 'AGC': 0,  
                    'AUA': 0, 'ACA': 0, 'AAA': 0, 'AGA': 0,  
                    'AUG': 0, 'ACG': 0, 'AAG': 0, 'AGG': 0,  
                    
                    'GUU': 0, 'GCU': 0, 'GAU': 0, 'GGU': 0,  
                    'GUC': 0, 'GCC': 0, 'GAC': 0, 'GGC': 0,  
                    'GUA': 0, 'GCA': 0, 'GAA': 0, 'GGA': 0,  
                    'GUG': 0, 'GCG': 0, 'GAG': 0, 'GGG': 0,  
                    }
                    
        self.dnaCodonCompDict = {key.replace('U','T'): value for key, value in rnaCodonCompDict.items()}

        
    def addSequence (self, inSeq):
        upperInSeq = inSeq.upper() #capitalized the string 

        for nucleotides in self.NucCompDict:
            if nucleotides in self.NucCompDict:
                self.NucCompDict[nucleotides] += 1
        for codon in range(0, len(codon), 3):
            if codon in self.dnaCodonTable.keys(): 
                self.dnaCodonCompDict[codon] += 1 #update in dnaCompDict
                self.aaCompDict[self.dnaCodonCompDict[codon]] += 1 #update in aaCompDict
            if codon in self.rnaCodonTable.key(): 
                self.rnaCodonCompDict[codon] += 1 #update in rnaCompDict
                self.aaCompDict[self.rnaCodonCompDict[codon]] += 1 #update in aaCompDict 

    def aaComposition(self):
        return self.aaCompDict
    def nucComposition(self):
        return self.nucCompDict
    def codonComposition(self):
        return self.codonCompDict
    def nucCount(self):
        return sum(self.nucComp)

class ProteinParam :
    '''
    Create class ProteinParam that includes methods that analyze the user input protein sequence.

    These tables are for calculating:
        aa2mw: molecular weight of each amino acid. 
        mwH2O: molecular weight of H2O.
        aa2abs280: absorbance at 280 nm.
        aa2chargePos: pKa of positively charged Amino Acids. 
        aa2chargeNeg: pKa of negatively charged Amino acids. 
        aaNterm: pKa of the N terminis.
        aaCterm: pKa of the C terminis.
    
    As written, these are accessed as class attributes.
    For example: ProteinParam.aa2mw['A'] or ProteinParam.mwH2O
    '''
    aa2mw = { 
            'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
            'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
            'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
            'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
            }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125} 

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6} 
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10} 
    
    aaNterm = 9.69 
    aaCterm = 2.34 

    def __init__ (self, protein):
        '''Creates a amino acid composition dictionary and computes the number of of each amino acid.'''
        upperProtein = protein.upper() #Upper cased sequence. 

        #Dictionary for every amino acid, starts with 0
        self.aaCompDict = {
                            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0, 
                            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
                            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
                            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0
                            }
        
        #For every character in upperProtein (input), if that character is in the dictionary, plus 1.
        #This counts how many valid characters there're in the sequence + disregards invalid characters. 
        for aa in upperProtein: 
            if aa in self.aaCompDict: 
                self.aaCompDict[aa] += 1 

    def aaCount (self):
        '''Computes and returns the length of the aa sequence by summing all values in aa dictionary.'''
        aaSeqLenth = sum(self.aaCompDict.values()) #all values in aa dictionary = total length of aa sequence. 
        return aaSeqLenth 
    
    def pI (self):
        '''Calculates and returns the theoretical isolelectric point. which is the pH that yields a neutral net charge.'''
        currentCharge = 999
        currentpH = 999
       
        for pH in range(0 ,1400): #iterates every pH in range 0-14
            pH = pH/100 
            newCharge = self._charge_(pH) #gets the charge value from _charge_ method
            if newCharge < currentCharge: #if charge used in the charge method is less than current charge,
                if newCharge > 0: #and above 0,
                    currentCharge = newCharge #then the value will be current charge
                    currentpH = pH #and the pH used in the charge method is the current pH 

        return currentpH

    def aaComposition (self) :
        '''Returns all values in the aa composition dictionary.'''
        return self.aaCompDict

    def _charge_ (self, pH):
        '''Calculates and returns the net charge on the protein at a specific pH.'''
        
        #Calculates net Positive charge. 
        netPosCharge = 0 #Inititalized at 0.
        for aa in ProteinParam.aa2chargePos.keys(): #adding 
            netPosCharge += self.aaCompDict.get(aa) * ((10 ** ProteinParam.aa2chargePos.get(aa))\
                          / (10 ** ProteinParam.aa2chargePos.get(aa) + 10 ** pH))
        netPosCharge += 10 ** ProteinParam.aaNterm / (10 ** ProteinParam.aaNterm + 10 ** pH)

        #Calculates net Negative charge.
        netNegCharge = 0
        for aa in ProteinParam.aa2chargeNeg.keys():
            netNegCharge += self.aaCompDict.get(aa) * ((10 ** pH)\
                          / (10 ** ProteinParam.aa2chargeNeg.get(aa) + 10 ** pH))
        netNegCharge += 10 ** pH / (10 ** ProteinParam.aaCterm + 10 ** pH)

        #Calculates and returns total net charge.
        netCharge = netPosCharge - netNegCharge
        return netCharge

    def molarExtinction (self): 
        '''Calculates and returns the molar extinction coefficient.'''   

        #molar extinction coefficient equation. the value of Y,W,C in aa dict times the abs value at 280nm.     
        molarExtinct = (self.aaCompDict.get('Y') * ProteinParam.aa2abs280.get('Y'))\
                     + (self.aaCompDict.get('W') * ProteinParam.aa2abs280.get('W'))\
                     + (self.aaCompDict.get('C') * ProteinParam.aa2abs280.get('C')) 
        return molarExtinct

    def massExtinction (self):
        '''Calculates and returns the Mass extinction coefficient.'''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0 

    def molecularWeight (self):
        '''Calculates and returns the molecular weight of the protein sequence.'''
        aaMolarWeight = 0 #mw initialized at 0 
        
        #for every amino acid in aa dictionary, the value of each aa is mutiplied with mw of the aa
        #minus the mw of water that are released with peptide bond formation.  
        for aa in self.aaCompDict.keys(): 
            aaMolarWeight += self.aaCompDict[aa] * (ProteinParam.aa2mw[aa] - ProteinParam.mwH2O)  
        
        finalMolarWeight = ProteinParam.mwH2O + aaMolarWeight
        return finalMolarWeight

import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
