#!/usr/bin/env python
# coding: utf-8

# # CLASS FILEREADER
# ---
# 
# This class reads the file passed in.
# 
# ### FUNCTIONS
# 
#     
#     init : class constructor

# In[4]:


class fileReader:

    def __init__(self, infile):
        '''
        Class constructor
        '''
        self.protein = ''
        self.data = ''

        with open(infile) as file:
            for line in file:
                self.data += line
                self.protein = line.rstrip()
        self.data = self.data.replace(self.protein, '').rstrip()
        self.pCodons = len(self.protein) * 3


# # CLASS AAFINDER
# ---
# 
# This class searches the genome and finds the matching polypeptides. It searches for a suitable beginning or ending codon, and then builds the rest of the string from there. 
# 
# ### FUNCTIONS AND ATTRIBUTES
# 
#     codons: a class attribute holding the codon table. 
#     init: class constructor
#     revComp: makes a reverse complement of a codon
#     findProteins: finds the polypeptide, if it exists
#     buildForward: builds the PP starting with the 5' starting AA
#     buildBackward: builds the PP starting from the 3' starting AA, working backward
#     printSequences: prints it all out
# 

# In[7]:





class AAfinder:
    
    #class attribute dictionary of codons
    codons = {'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'I': ['ATT', 'ATC', 'ATA'],
              'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
              'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
              'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'],
              'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
              'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG']}

    def __init__(self, infile):
        '''
        Class constructor
        '''
        self.fr = fileReader(infile)
        self.foundList = []



    # could be better done with a bunch of hardcoded str.replace in a line
    def revComp(self, codon):
        '''
        Function for reversing a codon
        
        Parameters
            codon - codon to reverse
        '''
        newString = ''
        for i in reversed(codon):
            if i == 'A':
                newString += 'T'
            elif i == 'T':
                newString += 'A'
            elif i == 'C':
                newString += 'G'
            elif i == 'G':
                newString += 'C'
        return newString



    def findProteins(self):
        '''
        Function that finds the polypeptide the in genome, if it exists
        '''
        #reads the genome to find a fitting starting point
        for x in range(len(self.fr.data) - len(self.fr.protein)):
            window = self.fr.data[x: x + 3]
            if window in AAfinder.codons[self.fr.protein[0]]:
                self.buildForward(x)
            elif self.revComp(window) in AAfinder.codons[self.fr.protein[-1]]:
                self.buildBackward(x)
            else:
                continue

    def buildForward(self, startIndex):
        '''
        Builds the PP from the 5' end
        
        Parameters
            startIndex - index in the genome to start at
        '''
        slices = []
        index = 0
        
        #sees if the rest of the PP can be made; if not, we move on
        for x in range(startIndex, startIndex + self.fr.pCodons, 3):
            if self.fr.data[x: x+3] in AAfinder.codons[self.fr.protein[index]]:
                slices.append(self.fr.data[x: x+3])
                index += 1
            else:
                slices.clear()
                break
        if len(slices) != 0:
            self.foundList.append(''.join(slices))


    def buildBackward(self, startIndex):
        '''
        Builds the PP from the 3' end instead
        
        Parameters
            startIndex - index in the genome to start at
        '''
        slices = []
        index = -1
        
        #sees if the rest of the PP can be made; if not, we move on
        for x in range(startIndex, startIndex + self.fr.pCodons, 3):
            if self.revComp(self.fr.data[x: x+3]) in AAfinder.codons[self.fr.protein[index]]:
                slices.append(self.fr.data[x: x+3])
                index -= 1
            else:
                slices.clear()
                break
        if len(slices) != 0:
            self.foundList.append(''.join(slices))


    def printSequences(self):
        '''
        Prints out all the sequences that were found to match the PP
        '''
        with open('p16out.txt', 'w') as out:
            for x in self.foundList:
                out.write(x + '\n')








# # MAIN
# This is the main function for instantiation and execution
# 

# In[8]:


def main(infile):
    aa = AAfinder(infile)
    aa.findProteins()
    aa.printSequences()



if __name__ == '__main__':
    main('rosalind_ba4b-3.txt')


# # INSPECTIONS
# ---
# 
# ### INSPECTOR: QUIN LAMOTHE
# 1. Concise readable code, efficient sequence parsing-- needs comments. **I FELT LIKE TOO MANY COMMENTS DECREASE READABILITY, ESPECIALLY IN A SHORT PROGRAM SUCH AS THIS**

# In[ ]:




