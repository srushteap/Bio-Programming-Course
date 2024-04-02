#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# Name: Srushti Patil (spatil5)
# Group Members: Dylan Joedicker, Andres Ricardez, Valerie Lentine
'''
Reads a file of fasta seqeunces from STDIN and then finds the unique subsequences
that occur in each single tRNA. Each of the 22 sets are minimized so that no subsequence set 
is a substring of any other part of the set. 
'''
import sys
import sequenceAnalysis

class findUnique:
    def __init__(self, tRNA_header, tRNA_sequences):
        '''Initialize with tRNA names and sequences.'''
        self.tRNA_header = tRNA_header
        self.tRNA_sequences = tRNA_sequences
        self.unique_list = [] #list to store unique substrings for each tRNA sequence
        self.essential_list = [] #list to store essential substrings for each tRNA sequence

    def powerSets(self):
        '''create powerset of sequences. power set is all the possible substrings of each length in a string'''
        self.powerset_list = [] # a list of all power sets 
        for sequence in self.tRNA_sequences:
            powerset = set()
            for index in range(len(sequence)):
                length = len(sequence)
                while length > index:
                    powerset.add(sequence[index : length]) #add all substrings starting at 'index'
                    length -= 1
            self.powerset_list.append(powerset)

    def uniques(self):
        '''Compute unique sets from power sets.'''
        for power_set in self.powerset_list:
            union_set = set() #set to store union of all other sequences' powersets
            for other_power_set in self.powerset_list:
                if other_power_set is not power_set:
                    union_set = union_set.union(other_power_set)
            unique_set = power_set.difference(union_set) #subtract to find unique substrings
            self.unique_list.append(unique_set)
        
    def essential_sets(self):
        '''compares subsets to eachother and keeps only the essential ones'''
        for unique_set in self.unique_list:
            non_essential = set() #set to store non-essential substrings
            for subString in unique_set:
                if subString[1:] in unique_set or subString[:-1] in unique_set:
                    non_essential.add(subString)
            essential_set = unique_set.difference(non_essential) #set containing all the essential substrings
            self.essential_list.append(essential_set)
        
    def printReport(self):
        Counter = 0 #tracks index of headers and sequences
        output = [] #list to store formatted output 
        
        for essential_set in self.essential_list:
            name = self.tRNA_header[Counter]
            sequence = self.tRNA_sequences[Counter]

            formattedessentials = [] #stores all essential substrings after being formatted
            for items in essential_set:
                itemLength = len(items)
                locations = [] #store start positions of substrings in the sequence
                for i in range(len(sequence) - itemLength + 1):
                    if sequence[i:i+itemLength] == items:
                        locations.append(i)
                
                for location in locations: 
                    dots = ('.') * location
                    formattedessentials.append(dots + items)

            formattedessentials.sort(key=lambda i: len(i)) #sorts formatted substrings by length to improve readability
            output.append([name, sequence, formattedessentials])
            Counter += 1
        
        for finaloutput in sorted(output): 
            print(finaloutput[0]) #print tRNA header
            print(finaloutput[1])#print tRNA sequence
            
            for output1 in finaloutput[2]: 
                print(output1) # Print each essential substring formatted with leading dots
            
def main():
    '''Takes sequences from stdin and prints report'''
    fasta_sequence = sequenceAnalysis.FastAreader()
    tRNA_header = []
    tRNA_sequences = []

    for headers, sequences in fasta_sequence.readFasta():
        tRNA_header.append(headers.replace("-", "").replace("_", "").replace(".", "").replace(" ", ""))
        tRNA_sequences.append(sequences.replace("-", "").replace("_", "").replace(".", "").replace(" ", ""))
     
    tRNAUnique = findUnique(tRNA_header, tRNA_sequences)
    tRNAUnique.powerSets()
    tRNAUnique.uniques()
    tRNAUnique.essential_sets()
    tRNAUnique.printReport()

if __name__ == "__main__":
    main()
