#!/usr/bin/env python3
# Name: Srushti Patil(spatil5)
# Group Members: Dylan Joedicker, Valerie Lentine, Andres Ricardez
'''
Calculates the physical-chemical properties of a user inputted protein sequence string. 
Example:
 Input: VLSPADKTNVKAAW
Output: Number of Amino Acids: 14
        Molecular Weight: 1499.7
        molar Extinction coefficient: 5500.00
        mass Extinction coefficient: 3.67
        Theoretical pI: 9.88
        Amino acid composition:
        A = 21.43%
        C = 0.00%
        D = 7.14%
        E = 0.00%
        F = 0.00%
        G = 0.00%
        H = 0.00%
        I = 0.00%
        K = 14.29%
        L = 7.14%
        M = 0.00%
        N = 7.14%
        P = 7.14%
        Q = 0.00%
        R = 0.00%
        S = 7.14%
        T = 7.14%
        V = 14.29%
        W = 7.14%
        Y = 0.00%
'''
class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    def __init__ (self, protein):
        '''Create a dictionary that holds aaComposition'''
        self.protein = protein
        self.aaDictionary = {aa: 0 for aa in 'ACDEFGHILKMNPQRSTVYW'} #sets each count in the dictionary as 0

    def aaCount (self):
        '''Counts the amino acids in the protein'''
        for aa in self.protein:
            aa = aa.upper() #capitalizes each amino acid to count for lower case inputted sequence
            if aa in self.aaDictionary:
                self.aaDictionary[aa] += 1
        return sum(self.aaDictionary.values())

    def pI (self, precision = 2):
        '''Finds pH value thats yields neutral net charge'''  
        pH_values = [i / 100 for i in range(1, 1401)] #generates a list of pH values from 0.01 to 14.00
        closest_pH = min(pH_values, key=lambda init_pH: abs(self._charge_(init_pH))) #finds the pH value that results in the closest value to zero
        return round(closest_pH, precision) #returns rounded pH value to specific precision 
        
    def aaComposition (self):
        '''Returns dictionary of amino acids with corresponding counts'''
        return self.aaDictionary

    def _charge_ (self, pH):
        '''Calculates net charge of protein at specific pH'''
        self.aa2chargePos = {'K': 10.5, 'R': 12.4, 'H':6}
        self.aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
        self.aaNterm = 9.69
        self.aaCterm = 2.34
        
        chargePos = 0
        #totals the positive charges
        for aa, i in self.aa2chargePos.items():
            chargePos += ( self.aaDictionary[aa] * (10 ** i) / ((10 ** i) + (10 ** pH)) )
        chargePos += ( (10 ** self.aaNterm) / ((10 ** self.aaNterm) + (10 ** pH)) )
        
        chargeNeg = 0
        #totals the negative charges
        for aa, i in self.aa2chargeNeg.items():
            chargeNeg += ( self.aaDictionary[aa] * (10 ** pH) / ((10 ** i) + (10 ** pH)) )
        chargeNeg += ( (10 ** pH) / ((10 ** self.aaCterm) + (10 ** pH)) )
       
        #finds the net charge
        netCharge = chargePos - chargeNeg
        
        return netCharge

    def molarExtinction (self):
        '''Calculates extinction coefficient'''
        self.aa2abs280 = {'Y':1490, 'W': 5500, 'C': 125}
        extinct = sum(self.aaDictionary[aa] * self.aa2abs280.get(aa,0) for aa in self.aa2abs280) #Use get method to handle cases where amino acid is not in aa2abs280
        return extinct
        
    def massExtinction (self):
        '''Calculates the mass extinction coefficient'''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''Calculates and returns molecular weight of the protein'''
        aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        } 
        mwH2O = 18.015
        
        mw = sum(aa2mw[aa] for aa in self.protein if aa in aa2mw)
        mw = mw - mwH2O * (len(self.protein)-1)
        return mw

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()