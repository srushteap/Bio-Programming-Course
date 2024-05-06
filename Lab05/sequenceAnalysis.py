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
        closest_pH = min(pH_values, key=lambda init_pH: abs(self._charge_(init_pH))) #finds the pH value that results in the closest           value to zero
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
        extinct = sum(self.aaDictionary[aa] * self.aa2abs280.get(aa,0) for aa in self.aa2abs280) #Use get method to handle cases where         amino acid is not in aa2abs280
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
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''inString (str): optional input sequence. defaults to an empty string
           initializes data structures to store counts for nucleotides, codons, and amino acids
        '''
        self.nucRNAComp = {'A': 0, 'C': 0, 'G': 0, 'U': 0, 'N': 0}
        self.nucDNAComp = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        self.codonComp = {} 
        self.aaComp = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0,
                       'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'U': 0,
                       'V': 0, 'W': 0, 'Y': 0, '-': 0}
        self.inString = inString.upper()
        
    def addSequence (self, inSeq):
        '''
        adds a new sequence to the existing data
        - in seq(str): input sequence of nucleotides {ACGTUN}
        
        This method must accept new sequences, from the {ACGTUN} alphabet, 
        and each sequence can be presumed to start in frame 1. This data must be added to the data that 
        you were given with the init method (if any).
        '''
        
        self.inString += inSeq.upper()
        return self.inString
    
    def aaComposition(self):
        ''' 
        Calculates and returns the amino acid composition
        
        Returns a dictionary with counts for each amino acid
        '''
        #iterate through codonComp and use translation table to update aaComp
        for i in range(0, len(0, len(self.inString) - 2, 3)):
            codon = self.inString[i:i+ 3] 
            
            if 'T' in codon:
                if codon in self.dnaCodonTable: 
                    self.aaComp[self.dnaCodonTable[codon]] += 1
            else:
                if codon in self.rnaCodonTable:
                    self.aaComp[self.rnaCodonTable[codon]] += 1 
       
        return self.aaComp
    
    def nucComposition(self):
        ''' 
        Returns the count of valid nucleotides. Distinguish between RNA and DNA nucleotides. N bases should also be counted
        '''  
        for i in self.inString:
            if i in 'ACGTUN':
                if 'U' in self.inString:
                    self.nucRNAComp[i] += 1
                elif 'T' in self.inString:
                    self.nucDNAComp[i] += 1              
                  
        return self.nucRNAComp, self.nucDNAComp
    
    def codonComposition(self):
        ''' 
        Returns counts of codons. Codons with invalid bases should be discarded. codons with N should be discarded. 
        All counts stored as RNA. 
        '''
        newSeq = self.inString.replace('T', 'U')
            
        for i in range(0,len(newSeq) - 2, 3):
            codon = self.inString[i:i + 3]
            if 'N' in codon:
                continue
            if all (base in "ACGTUN" for base in codon):
                self.codonComp[codon] = self.codonComp.get(codon, 0) + 1       
        
        return self.codonComp
    
    def nucCount(self):
        ''' 
        returns count of every valid nucleotide
        '''
        nucSum = self.nucRNAComp.values() + self.nucDNAComp.values()
        return nucSum

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
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
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

class OrfFinder:
    '''
    Takes in a sequence of a FASTA file and stores ORFs in a lists of dictionaries.
    '''
    # Define start and stop codons and the complement dictionary for DNA sequences
    start_codons = ["ATG"]
    stop_codons = ["TAG", "TAA", "TGA"]
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, seq):
        # Initialize OrfFinder object with the DNA sequence
        self.seq = seq
        self.orfs = []  # List to store identified ORFs

    def findOrfs(self):
        '''Find open reading frames (ORFs) in the provided DNA sequence.'''
        start_position = []  # Keeps track of the start codon positions
        foundStart = 0
        foundCodon = 0

        for frame in range(0, 3):  # Iterate through three reading frames
            foundStart = 0
            foundCodon = 0
            start_position = []  # Clear start position list for each frame
            for i in range(frame, len(self.seq), 3):  # Iterate over the sequence in sets of three nucleotides
                codon = self.seq[i:i + 3]  # Codon length is 3 nucleotides
                if codon == 'ATG':
                    start_position.append(i)
                    foundStart = 1
                    foundCodon = 1
                if codon in OrfFinder.stop_codons and foundStart:
                    start = start_position[0] + 1 - frame
                    stop = i + 3
                    length = stop - start + 1
                    self.saveOrf((frame % 3) + 1, start, stop, length)
                    start_position = []
                    foundStart = 0
                    foundCodon = 1
                if not foundCodon and codon in OrfFinder.stop_codons and foundStart:
                    start = 1
                    stop = i + 3
                    length = stop - start + 1
                    self.saveOrf((frame % 3) + 1, start, stop, length)
                    start_position = []
                    foundCodon = 1
            if foundStart:
                start = start_position[0] + 1
                stop = len(self.seq)
                length = stop - start + 1
                self.saveOrf((frame % 3) + 1, start, stop, length)
        return self.orfs

    def revComp(self):
        '''Return the reverse complement of the DNA sequence.'''
        revSeq = ''.join([self.complement[base] for base in self.seq[::-1]])
        return revSeq

    def findRevOrfs(self):
        '''Find open reading frames (ORFs) in the reverse complement of the DNA sequence.'''
        reverse = self.revComp()  # Get the reverse complement
        start_position = []
        foundStart = 0
        foundCodon = 0

        for frame in range(0, 3):  # Iterate through three reading frames
            foundStart = 0
            foundCodon = 0
            start_position = []  # Clear start position list for each frame
            for i in range(frame, len(reverse), 3):  # Iterate over the sequence in sets of three nucleotides
                codon = reverse[i:i + 3]  # Codon length is 3 nucleotides
                if codon == 'ATG':
                    start_position.append(i)
                    foundStart = 1
                    foundCodon = 1
                if codon in OrfFinder.stop_codons and foundStart:
                    stop = len(reverse) - start_position[0]
                    start = len(reverse) - (i + 2)
                    if frame == 1:
                        stop += 1
                    elif frame == 2:
                        stop += 2
                    length = stop - start + 1
                    self.saveOrf(-1 * ((frame % 3) + 1), start, stop, length)
                    start_position = []
                    foundStart = 0
                    foundCodon = 1
                if not foundCodon and codon in OrfFinder.stop_codons and foundStart:
                    start = len(reverse) - i - 2
                    stop = len(reverse)
                    length = stop - start + 1
                    self.saveOrf(-1 * ((frame % 3) + 1), start, stop, length)
                    start_position = []
                    foundCodon = 1
            if foundStart:
                start = start_position[0] + 1
                stop = 1
                length = stop - start + 1
                self.saveOrf(-1 * ((frame % 3) + 1), start, stop, length)
        return self.orfs

    def saveOrf(self, frame, start, stop, length):
        """Save ORF information."""
        self.orfs.append({'frame': frame, 'start': start, 'stop': stop, 'length': length})
