'''
Reads sequence in a file and returns sequence length, GC content, and amino acid content. 
'''
def main (fileName=None):
    myReader = FastAreader(fileName) 
    myNuc = NucParams()
    for head, seq in myReader.readFasta(): #read file or stdin and add sequences from file to NucParams
        myNuc.addSequence(seq)
        
    '''sort codons in alphabhetical order, by Amino Acid'''
    seq_length_mb = myNuc.nucCount() / 1e6 #calculate sequence length in mb 
    print(f'sequence length = {seq_length_mb:.2f} Mb\n') 
    
    #calculate and print GC content
    gc_content = (myNuc.nucComposition()["G"] + myNuc.nucComposition()["C"]) / myNuc.nucCount() * 100
    print(f'GC content = {gc_content:.1f}%\n')
    
    #sorts codons alphabetically by amino acid
    nucs = sorted(myNuc.codonComposition().keys(), key=lambda x: NucParams.rnaCodonTable[x])
    
    #calculates and prints relative codon usage
    for nuc in nucs:
        aa = NucParams.rnaCodonTable[nuc]
        val = myNuc.codonComposition()[nuc] / myNuc.aaComposition()[aa]
        this_codon_comp = myNuc.codonComposition()
        print('{:s} : {:s} {:5.1f} ({:6d})'.format(nuc, aa, val * 100, this_codon_comp[nuc]))

if __name__ == "__main__":
    main() 
    