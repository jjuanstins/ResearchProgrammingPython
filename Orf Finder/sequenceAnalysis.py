#!/usr/bin/env python3
# Name: Justin Juang (jjuang)
# Group Members: None

'''
Program sequenceAnalysis.py is essentially a toolkit that inclues modules that are helpful
in analysis of biological sequences, including DNA, RNA, and protein sequences. 

Modules include:
    OrfFinder: finds and stores open reading frames in a sequence.
    NucParams: the module counts the number of nucleotides, codons, and amino acids of a given 
        sequence and stores the information in respective dictionaries. 
    ProteinParams: the module analyzes a given protein sequence, and returns the information.
    FastAReader: the module opens and reads a FastA file.

Personal Note:
Lab 3: ProteinParams
Lab 4: NucParams + FastAReader
Lab 5: OrfFinder
'''

################################################################################################################################################

class OrfFinder():
    '''Creates class OrfFinder that finds ORF in a given sequence.'''

    def __init__(self, seq):
        '''
        Init method that creates objects for later use.
            orfList is the list of orf found.
            complementStrand is the complement strand sequence. 
        '''
        
        self.seq = seq #no need for upper() since file came with capitalized sequence

        self.orfList = [] #lists to store found orfs
        
        self.startCodons = ['ATG'] #list of start codon (no extra credit)
        self.stopCodons = ['TGA', 'TAG', 'TAA'] #list of stop codons
        
        #personal note: str.maketrans and .translate() goes together. 
        nucComplements = str.maketrans('ATGC', 'TACG') #cool function i found str.maketrans on google that makes a translation table
        self.complementStrand = self.seq.translate(nucComplements)[::-1] #.translate is the function that actually translates

    def forwardFrame(self):
        '''Finds ORFs in the forward frame.'''

        startList = [] #list of all found start codon positions
        startFound = 0 #indicates if start codon is found 
        codonFound = 0 #indicates if a codon is found 

        for frame in range(0, 3): #for each frame (3 nucleotides at a time)  
            
            startFound = 0
            codonFound = 0 #indicators initialized at 0
            startList = [] #start codon list cleared
            
            for nucleotide in range(frame, len(self.seq), 3): #for each nucleotide (poisition) in that frame
                
                codon = self.seq[nucleotide : nucleotide + 3] #a codon is three nucleotides
                
                #if the codon is in start codon (or ATG)
                if codon in self.startCodons: 

                    startList.append(nucleotide) #adds the position of the nucleotide to the start list ("A" pos)
                    startFound = 1 #indicates that a start codon is found 
                    codonFound = 1

                #if codon is in stop codons and startFound = 1 (normal orf)
                if codon in self.stopCodons and startFound == 1:

                    startPos = startList[0] + 1 - frame #the start pos of this orf is the first position in the list  
                    stopPos = nucleotide + 3 #and stop pos is at the nucleotide position + 3 (for codon) 
                    orfLength = stopPos - startPos + 1 #and orf length is just pos of stop - start
                    
                    orf = [(frame%3) + 1, startPos, stopPos, orfLength] #define orf
                    self.orfList.append(orf) #and append it 
                    
                    startList = [] #clears start list 
                    startFound = 0 
                    codonFound = 1 

                #dangling stop: no start codon found (lonely stop dangling there without start)
                if codonFound == 0 and codon in self.stopCodons: 
                    startPos = 1 #the start pos would be the start of sequence (so 1)
                    stopPos = nucleotide + 3 #and end pos would be the pos of stop codon + 3 
                    orfLength = stopPos - startPos + 1 #orf length
                    
                    orf = [(frame%3) + 1, startPos, stopPos, orfLength] #define orf
                    self.orfList.append(orf) #and append it

                    startList = [] #clear start list 
                    codonFound = 1

            #dangling start: no stop found (but start found)
            if startFound == 1:
                startPos = startList[0] + 1 #then start pos is the pos of start codon
                stopPos = len(self.seq) #and stop pos is the end of sequence (no stop)
                orfLength = stopPos - startPos + 1
                orf = [(frame%3) + 1, startPos, stopPos, orfLength]
                self.orfList.append(orf)

        return self.orfList #returns orf list with found orfs 

    def reverseFrame(self): 
        '''Finds ORFs in the reverse frame on the complement frame.'''  
        
        #bascially same as forward frame... with minor differences
        startList = [] #empty start codon list
        startFound = 0
        codonFound = 0

        for frame in range(0, 3): #iterates through each frame (3 nucleotides)
            startFound = 0  #all initializes at 0
            codonFound = 0
            startList = []  

            for nucleotide in range(frame, len(self.complementStrand), 3): #iterates first position at each frame
                
                codon = self.complementStrand[nucleotide : nucleotide + 3] #creates codon

                #if start codon is found 
                if codon in self.startCodons:  
                    
                    startList.append(nucleotide) #append the start codon position
                    startFound = 1
                    codonFound = 1

                #if stop codon is found and there's a start codon
                if codon in self.stopCodons and startFound == 1: 
                    
                    stopPos = len(self.complementStrand) - startList[0] #gets stop position
                    startPos = len(self.complementStrand) - (nucleotide + 2) #gets start position
                    if frame == 1: 
                        stopPos += 1
                    elif frame == 2: 
                        stopPos += 2
                    
                    orfLength = stopPos - startPos + 1 #gets orf length
                    orf = [-1 * ((frame%3) + 1), startPos, stopPos, orfLength] #creates orf
                    self.orfList.append(orf) #append to orf list 

                    startList = [] #resets start list
                    startFound = 0
                    codonFound = 1

                #dangling stop: no start codon found 
                if codonFound == 0 and codon in self.stopCodons: 
                    
                    startPos = len(self.complementStrand) - nucleotide - 2 #gets start length
                    stopPos = len(self.complementStrand) #gets stop length (length of complement strand)
                    
                    orfLength = stopPos - startPos + 1 #gets orf length
                    orf = [-1 * ((frame%3) + 1), startPos, stopPos, orfLength] #creates orf 
                    self.orfList.append(orf) #appends orf

                    startList = []
                    codonFound = 1

            #dangling start: no stop found (but start found)
            if startFound == 1: 
                startPos =  startList[0] + 1 #get start pos 
                stopPos = 1 #stop is just the beginning of sequence
                orfLength = stopPos - startPos + 1 #gets length
                orf = [-1 * ((frame%3) + 1), startPos, stopPos, orfLength] #creates orf
                self.orfList.append(orf) #append orf 

        return self.orfList

################################################################################################################################################

class NucParams:
    '''Create class NucParams stores the number of amino acid, codon, and nucleotides in dictionaries '''
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
                    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
                    }
    dnaCodonTable = {key.replace('U','T'): value for key, value in rnaCodonTable.items()} #dna codon table

    def __init__ (self, inString=''):
        '''
        Instantiates the amino acid composition, nucleotide composition, and RNA codon composition 
        dictionaries. 
        '''
        self.aaCompDict = { #amino acid composition dictionary 
                            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0, 
                            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
                            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
                            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0,
                            '-': 0
                          }
       
        self.nucCompDict = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'U': 0, 'N': 0} #nucleotide composition dictionary

        self.rnaCodonCompDict = { #RNA codon composition dictionary
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

        self.dnaCodonCompDict = {key.replace('U','T'): value for key, value in self.rnaCodonCompDict.items()}
    
    def addSequence (self, inSeq):
        '''Accepts new sequences and stores nucleotide and RNA codon info in respective dictionaries.'''
        
        upperInSeq = inSeq.upper() #capitalizes the input string.

        #counts the nucleotides.
        for nucleotides in upperInSeq: 
            if nucleotides in self.nucCompDict:
                self.nucCompDict[nucleotides] += 1

        #counts RNA codons.
        for nucleotide in range(0, len(upperInSeq), 3): 
            codon = upperInSeq[nucleotide: nucleotide + 3] #defining codons
            rnaCodon = codon.replace('T','U')
            if rnaCodon in self.rnaCodonTable.keys(): 
                self.rnaCodonCompDict[rnaCodon] += 1 #update in rnaCompDict
                self.aaCompDict[self.rnaCodonTable[rnaCodon]] += 1 #update in aaCompDict 

    def aaComposition(self):
        '''Returns the amino acid composition dictionary.'''
        return self.aaCompDict
    def nucComposition(self):
        '''Returns the nucleotide composition dictionary.'''
        return self.nucCompDict
    def codonComposition(self):
        '''Returns the codon composition dictionary.'''
        return self.rnaCodonCompDict
    def nucCount(self):
        '''Returns the sum of all values in the nucleotide composition dictionary (Length of sequence).'''
        return sum(self.nucCompDict.values())


################################################################################################################################################

class ProteinParam :
    '''Create class ProteinParam that includes methods that analyze the user input protein sequence.'''
    
    aa2mw = { #molecular weight of each amino acid. 
            'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
            'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
            'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
            'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
            }

    mwH2O = 18.015 #molecular weight of H2O
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125} #absorbance at 280 nm.

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6} #pKa of positively charged Amino Acids
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10} #pKa of negatively charged Amino acids.
    
    aaNterm = 9.69 #pKa of the N terminis.
    aaCterm = 2.34 #pKa of the C terminis.

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

################################################################################################################################################

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
                if not line: # we are at EOF
                    return header, sequence
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

