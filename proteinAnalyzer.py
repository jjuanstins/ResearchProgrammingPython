#!/usr/bin/env python3
# Name: Justin Juang (jjuang)
# Group Members: None

from sequenceAnalysis import ProteinParam
import sys

def main():
    '''main method that prints out returned values'''
    inString = input('protein sequence? ')
    
    while inString:
        myParamMaker = ProteinParam(inString) #instance 
        myAAnumber = myParamMaker.aaCount() #get aa length
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition: ")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()