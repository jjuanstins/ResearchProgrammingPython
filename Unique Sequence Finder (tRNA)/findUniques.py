#!/usr/bin/env python3
# Name: Justin Juang (jjuang)
# Group Members: None - With assistance from tutor and TA.

'''
Program findOrfs.py reads a file of fasta sequences from STDIN, finds the unique subsequences
that occur in a single tRNA such that no members of this set occur among any of the other tRNA sets.

Program includes methods find uniquetRNA, essentialtRNA, fastAreader, and main:
    uniquetRNA: finds the unique tRNA sequences.
    essentialtRNA: finds the essential tRNA sequences from the unique sequences.
    fastAreader: reads a file of fasta sequences from STDIN.
    main: prints out header, sequence, and essential sequences of each tRNA     
'''


class FindUniques:
    '''Create class FindUnique that finds unique sequence in mitochondrial tRNA'''

    tRNAList = list() #create empty list to store tRNA elements

    def __init__(self, header, sequence):
        '''Creates power set for every tRNA sequence.'''

        self.sequence = sequence  
        self.header = header 

        self.powerSet = set() #creates empty set to store every element of power set from sequence 

        #creating power set
        for xNuc in range(len(sequence)): #for every nucleotide in the pos x sequence 
            for yNuc in range(xNuc, len(sequence)): #and from x to end (pos yNucleotide)
                self.powerSet.add(sequence[xNuc : yNuc]) #add sequence from x to y to powerset set (.add wont add if element already exist in set)

        FindUniques.tRNAList.append(self.powerSet) #appends power set to overall tRNA set 


    def uniquetRNA(self):
        '''Creates method uniquetRNA that finds and returns a set of all unique tRNA.'''

        thistRNA = set() #create set thistRNA
        uniquetRNAs = set() #create set uniquetRNA storing unique tRNAs

        for powerSet in FindUniques.tRNAList: #for every tRNA element in tRNASet (power set)
            if powerSet is not self.powerSet:
                thistRNA.update(powerSet)

            uniquetRNAs = self.powerSet - thistRNA #taking the difference will yield unique sequences. 

        return uniquetRNAs #return unique tRNAs 

    def essentialtRNA(self): 
        '''Creates method essentialtRNA that finds and returns a set of all essential tRNA'''

        nonEssentials = set() #create empty set for all non essential tRNAs 
        essentialtRNAs = set() #create empty set for all eseential tRNAs

        uniqueSet = self.uniquetRNA()
        sortedUniques = sorted(list(uniqueSet), key=len, reverse = False) #sort unique tRNAs by len, shortest first

        #finds essential - adds non essential to set and take the difference later.
        for uniques in sortedUniques: 
            for otherUniques in sortedUniques: 
                if otherUniques is not uniques: 
                    if uniques in otherUniques: 
                        nonEssentials.add(otherUniques) 

        essentialtRNAs = uniqueSet - nonEssentials #difference

        return essentialtRNAs #returns essential tRNAs (set)

import sys
class FastAreader:
    '''Create method FastAreader that reads a fasta file from STDIN.'''
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        '''Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        '''Read an entire FastA record and return the sequence header/sequence.'''
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

#############################################################################################
# Main
# Here is the main program
# 
#############################################################################################

def main(inCL=None):
    '''Method main that prints out header, sequence, and all uniques/essential sequences for each tRNA sequence.'''
   
    myFastAReader = FastAreader() #instantiating FastAReader class
    tRNAList = list() #create tRNA List 
    
    for header, sequence in myFastAReader.readFasta(): #for every header and sequence that are parsed from readFastA method
        myFindUniques = FindUniques(header, sequence) #instantiating FindUniques class
        tRNAList.append(myFindUniques)
    
    #sortedtRNAList = sorted(sortableHeader)

    headerSortedDotList = [] #[[header][sortedDot]]

    for myFindUniques in tRNAList: 
        #sortedtRNAList = sorted(tRNAList)

        #print(myFindUniques.header) #prints header
        #print(myFindUniques.sequence) #prints original sequence
        
        myFindUniques.uniquetRNA() #instantiating uniquetRNA method
        myFindUniques.essentialtRNA() #instantiating essentialtRNA method

        #calculating/obtaining dot numbers before essential sequences 
        dotSet = list() 

        for tRNA in myFindUniques.essentialtRNA(): 
            dotPos = myFindUniques.sequence.find(tRNA)
            dotSet.append(("." * dotPos) + tRNA)

        sortedNumDot = sorted(dotSet, key = len, reverse = False)
        # print(sortedNumDot)
        headerSortedDotList.append([myFindUniques.header,sortedNumDot])
        # print(headerSortedDotList)

    finalOutput = sorted(headerSortedDotList) 
    #prints sequence
    for item in finalOutput:
        # print(item)     
        print(item[0])
        for seq in item[1]:
            print(seq)

if __name__ == "__main__":
    main()  
