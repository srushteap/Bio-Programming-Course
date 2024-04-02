#!/usr/bin/env python3 
# Name: Srushti Patil(spatil5)  
# Group Members: Andres, Valerie, Dylan

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str): #DNAstring is defined as a kind of python string
  
  def length (self): #length method returns length of string
    return (length(self))
  def purify(self): 
    ''' Returns an upcased version of the string, collapsing a single run of Ns.'''
    st = '' 
    count = 0
    for i in range(0,len(self)): #iterates through the string and counts the number of Ns assuming that the Ns are together in one block
        if self[i] == 'N':
            st += 'N'   #string containing block of Ns
            count += 1
    self = self.replace(st,("{" +str(count) + "}")) #replaces block of Ns with count of Ns
    pureDNA = self.upper() #capitalizes string
    return pureDNA
      
def main(): 
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring(data) #creates a new name of class DNAstring 
    pureData = thisDNA.purify() #sends purify() method message to purify inputted DNA string
    print(pureData) #prints resulting purified DNA
     
main()