#!/usr/bin/env python3 
# Name: Srushti Patil(spatil5) 
# Group Members: Andres, Valerie, Dylan

'''
Collects a string of a codon or amino acid from user input. Looks up the information in the dictionaries, 
and returns the corresponding codon base or amino acid. 

Example:
 Input: ATG
Output: MET

 Input: UAG
Output: ---

 Input: E
Output: GLU

 Input: ASP
Output: D

'''
short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            }

long_AA = {value:key for key,value in short_AA.items()}

RNA_codon_table = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}
dnaCodonTable = {key.replace('U','T'):value for key, value in RNA_codon_table.items()}

def main():
    '''Asks user for codon or amino acid and searches through the codon and amino acids tables for corresponding value'''
    codonAA = input("Codon?") 
    codonAA = codonAA.upper() #not case sensitive, capitalizes input first
    
    exists = False #to check if input is valid
    
    #look for key and value for RNA and DNA codon tables
    for key in RNA_codon_table:
        if key == codonAA:
            print(codonAA + " = " + RNA_codon_table[key].upper())
            exists = True
        elif RNA_codon_table[key] == codonAA:
            print(RNA_codon_table[key].upper() + " = " + codon)
            exists = True
    for key in dnaCodonTable:
        if key == codonAA:
            print(codonAA + " = " + dnaCodonTable[key].upper())
            exists = True
        elif dnaCodonTable[key] == codonAA:
            print(dnaCodonTable[key].upper() + " = " + codon)
            exists = True
    
    #looks for key and value for amino acids
    for key in short_AA:
        if key == codonAA:
            print(codonAA + " = " + short_AA[key])
            exists = True
    for key in long_AA:
        if key == codonAA:
            print(codonAA + " = " + long_AA[key])
            exists = True
   
    #if not valid input, outputs "unknown"
    if exists != True:
        print("unknown")

main()