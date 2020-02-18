#!/usr/bin/env python3
# Name: Justin Juang (jjuang)
# Group Members: None

'''
Program findOrfs.py includes modules CommandLine and main.

CommandLine class creates a comand line program that allows user interface that is accessed 
through the terminal winder. 
Main class instantiates the OrfFinder module from sequenceAnalysis.py and prints out ORF info.

Required pseudocode: 
same reasoning with forward and backward frame
    for each frame(0-3)
        for each nucleotide from 0-length of sequence by 3
            define codon 
            if codon is a start codon
                then append that position
            if codon is a stop codon
                and if there's a previous start codon found
                    then make orf 
                if no start (dangling stop)
                    make orf from beginning of sequence
            if codon is a start but no stop is found (dangling start)
                then orf is from start codon position to end of sequence

'''

########################################################################
# CommandLine
########################################################################

from sequenceAnalysis import OrfFinder, FastAreader #imports OrfFinder and FastAReader class

class CommandLine():
    """
    Creates class CommandLine that allows user option input. 

    Description from Lab 5 pdf: 
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.
    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    """
    def __init__(self, inOpts=None):
        """
        Implement a parser to interpret the command line argv string using argparse.
            --> Argparse adjusted so that default is minGene = 100
        """
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Finds the largest ORF in a DNA sequence',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store_true', help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(100, 200, 300, 500, 1000), action='store',
                                 help='minimum Gene length', default = 100)
        self.parser.add_argument('-s', '--start', action='append', nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


########################################################################
# Main
# Here is the main program
#
#
########################################################################

def main(inCL=None):
    ''' 
    Creates function main that reads the FastA file, instantiates OrfFinder module and prints out 
    ORF info with user given options. 
    '''

    if inCL is None:
        myCommandLine = CommandLine()
        
        if myCommandLine.args.longestGene:
            fastaFile = FastAreader()
            
            for header, sequence in fastaFile.readFasta(): #reads FastA file
                print(header) #prints header
                
                myOrfFinder = OrfFinder(sequence) #institating class OrfFinder
                myOrfFinder.forwardFrame() #instantiating method forwardFrame 
                myOrfFinder.reverseFrame() #instantiating method reverseFrame
                
                #filters orf, only keeping orfs longer than minGene (user specified length)
                #personal note: turns out there's a lambda filter
                filteredOrf = filter(lambda orf:orf[3] > myCommandLine.args.minGene, myOrfFinder.orfList)
                
                #sorting orfs by decreasing length (lambda sort)
                #personal note: reverse = True reverses its. 
                sortedOrf = sorted(filteredOrf, key=lambda orf:orf[3], reverse = True) 
                
                #prints it out 
                for frame, startPos, stopPos, orfLength in sortedOrf:
                    print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, startPos, stopPos, orfLength))
    else:
        myCommandLine = CommandLine(inCL)
    print(myCommandLine.args)

if __name__ == "__main__":
    main()

