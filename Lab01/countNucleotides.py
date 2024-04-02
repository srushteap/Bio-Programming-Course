#!/usr/bin/env python3
# Name: Srushti Patil(spatil5)
# Group Members: none
'''
Builds a dnaString object and returns the length and number of bases of a DNA sequence.

input: a string of arbitrary length containing a DNA sequence 
output: the length of the DNA sequence and the count of the bases printed on the screen
'''
class dnaString (str): 
    def length (self): #length method returns length of string object from dnaString class
        return (len(self))
    
    def countNucleotideA (self): #creates method that counts number of As in object from dnaString class
        num_A = self.count('A')
        return num_A
    
    def countNucleotideC (self): #creates method that counts number of Cs in object from dnaString class
        num_C = self.count('C')
        return num_C
    
    def countNucleotideG (self): #creates method that counts number of Gs in object from dnaString class
        num_G = self.count('G')
        return num_G
    
    def countNucleotideT (self): #creates method that counts number of Ts in object from dnaString class
        num_T = self.count('T')
        return num_T

dna = input("Enter a dna sequence: ") #asks user to input DNA sequence
upperDNA = dna.upper() #capitalizes string
coolString = dnaString(upperDNA) #creates new name of class dnaString with capitalized string 

#prints length of DNA sequence and count of each base
print ("Your sequence is {0} nucleotides long with the following breakdown of bases:".format(coolString.length()))
print ("number As = {0} number Cs = {1} number Gs = {2} number Ts = {3}".format(
    coolString.countNucleotideA(),
    coolString.countNucleotideC(), 
    coolString.countNucleotideG(), 
    coolString.countNucleotideT() ) )
