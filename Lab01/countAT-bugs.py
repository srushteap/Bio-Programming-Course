#!/usr/bin/env python3
# Name: Srushti Patil(spatil5)
# Group Members: none
'''
Builds a dnaString object and returns the AT content of a DNA sequence. There were 3 bugs in the original code at lines 15, 17, and 23.

input: a string of arbitrary length containing a DNA sequence 
output: the AT content of the sequence printed on the screen
'''
class dnaString (str): 
    def length (self): #creates method that returns length of string object from dnaString class
        return (len(self))

    def getAT (self): #creates method that counts amount of AT in object from dnaString class
        num_A = self.count('A') #was not a string
        num_T = self.count("T")
        return ((num_A + num_T)/ len(self) ) #was wrong format

dna = input("Enter a dna sequence: ") #asks user to input dna sequence
upperDNA = dna.upper() 
coolString = dnaString(upperDNA) #creates new name of class dnaString with user input as its object
#prints AT content by sending message to getAT() to return AT content of the inputted string
print ("AT content = {0:0.3f}".format(coolString.getAT()) ) #was only 1 sig fig