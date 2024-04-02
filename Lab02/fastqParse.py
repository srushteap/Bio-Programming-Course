#!/usr/bin/env python3 
# Name: Srushti Patil(spatil5)
# Group Members: Andres, Valerie, Dylan

'''
Reads the seqname line of a FASTQ file from user input and parses out each field, returning a string for each field. 

Example:
 input: @EAS139:136:FC706VJ:2:2104:15343:197393
output: Instrument = EAS139
        Run ID = 136
        Flow Cell ID = FC706VJ
        Flow Cell Lane = 2
        Tile Number = 2104
        X-coord = 15343
        Y-coord = 197393
'''

class FastqString (str): #FastqString is now defined as a kind of python string
    ''' Class docstring goes here.'''
    def parse(self):
        '''Returns parsed out version of string. Displays each field with its corresponding value on a new line.'''
        self = self[1:] #first character is "@", removes "@" from string
        fileList = self.split(":") #each field is seperated by ":", splits string by ":" and creates a list holding each field
        
        '''Inputted string has 7 fields. Prints out each field and it's corresponding value'''
        print("Instrument = {0}".format(fileList[0]))
        print("Run ID = {0}".format(fileList[1]))
        print("Flow Cell ID = {0}".format(fileList[2]))
        print("Flow Cell Lane = {0}".format(fileList[3]))
        print("Tile Number = {0}".format(fileList[4]))
        print("X-coord = {0}".format(fileList[5]))
        print("Y-coord = {0}".format(fileList[6]))
    
def main():
    '''Gets user FASTQ file and parse its information'''
    file = input('FASTQ file seqname line?') #asks user for seqname line of FASTQ file 
    thisFile = FastqString(file) #creates a new name of FastqString class
    org = thisFile.parse() #sends parse() method to parse inputted file string

main()