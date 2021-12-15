#!/usr/bin/env python
# coding: utf-8

# # STRICTCYCLONE
# ---
# This class creates strictCyclone objects that search for possible peptides that match a given ideal spectrum. The algorithm used is the Branch and Bound algorithm, where all peptides are synthesized through expansion of all amino acids, then are excised if not meeting the desired criteria. 
# 
# ### FUNCTIONS/ATTRIBUTES
# 
# 
#     aminos - class attribute containing integer masses of amino acids
#     init - class constructor
#     lavarBurton - reads the file being used
#     extendoClip - does the branching and some bounding
#     checkSpecs - does rest of the bounding

# In[1]:



class strictCyclone:
    '''
    Class for finding peptides matching a given ideal spectrum
    '''
    
    aminos = {57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186}

    def __init__(self, infile):
        '''
        Class constructor
        
        Parameter
            infile - file to read in data from
        '''
        self.data = []  # set of ints that have fragment masses
        self.singles = []  # list single AA's as ints
        self.lavarBurton(infile)
        self.finalPeps = []
        self.spectrums = {}

    
    def lavarBurton(self, infile):
        '''
        Function to read the file and get the spectrum
        
        Parameter
            infile - file to read
        '''
        with open(infile) as file:
            
            for line in file:
                holdList = line.split()
                
                for weight in holdList:
                    self.data.append(int(weight))
                    if int(weight) in strictCyclone.aminos:
                        self.singles.append(int(weight))

    
    def extendoClip(self):
        '''
        Functions that extend the chain of amino acids
        '''
        oldPeps = [[i] for i in self.singles]
        newPeps = []
        
        #keeps extending the amino acid chains while there are still chains to extend
        while len(oldPeps) != 0:  
            
            for pep in oldPeps:  
                
                for nextA in self.singles:
                    if (pep[-1] + nextA) in self.data and sum(pep) < max(self.data):
                        newPeps.append([i for i in pep] + [nextA])
            oldPeps.clear()
            
            # checks if the fragments exist in the spectrum and adds them if so
            for newp in newPeps:
                if sum(newp) in self.data:
                    oldPeps.append(newp)
            newPeps = oldPeps
            
            # accounts for duplicates and fully formed chains
            for op in oldPeps:
                if sum(op) == max(self.data) and op not in self.finalPeps:
                    self.finalPeps.append(op)
                    newPeps.remove(op)
            oldPeps = newPeps
            newPeps = []

    #  checks if the spectrums of the synthesized peptides matches the parent spectrum
    def checkSpecs(self):

        window = 1
        holdList = []
        spectrum = [0]
        
        #creates a spectrum for each peptide found
        for pep in self.finalPeps:

            while window <= len(pep):
                for i in range(len(pep)):
                    for e in range(i, i + window):
                        try:
                            holdList.append(str(pep[e]))
                        except IndexError:
                            holdList.append(str(pep[(e - len(pep))]))
                    spectrum.append(sum([int(x) for x in holdList]))
                    if window == len(pep):
                        break
                    holdList.clear()
                window += 1
                
            #if the spectrums match, then add them to spectrums
            if sorted(spectrum) == sorted(self.data):
                self.spectrums['-'.join([str(i) for i in pep])] = spectrum
            window = 1
            spectrum = [0]
            holdList.clear()

        with open('p18out.txt', 'w') as out:
            for key in self.spectrums:
                out.write(key + '\n')




# # MAIN
# ---
# This is the main function for instantiation and execution
# 
# ### FUNCTIONS
# 
# 
#     main - function for executing the program
#     
#     

# In[4]:


def main(infile):
    sc = strictCyclone(infile)
    sc.extendoClip()
    sc.checkSpecs()
    
if __name__ == '__main__':
    main('rosalind_ba4e-2.txt')


# # INSPECTIONS
# ---
# ### Inspected By: Chris Condon
# 1. I like your naming conventions, and documentation is fairly easy to read
# 2. I think some bound checks, like checking for duplicates, may be redundant and slow your code down a bit. **WILL FIX LATER**
# 
# ### Inspected By: Aishwarya Basude
# 1. Interesting that you decided to do file reading in your class, i find it a little strange though because its not always the same infile format.
# 2. You don't seem to use the pseudocode as well so I'm not sure how efficient it is.
# 
# ### Inspected By: Drew Thompson
# 1. Curious why you named your file readin lavarBurton 
# 2. I find the structure of the code a little unintuitive for me personally (ie checking spectrums outside of your branch and bound), but it doesn't seem wrong
# 3. Comments are good and make it easier to follow
# 
# ### Inspected By: Immaad Mir
# 1. I like how you gave credit to Lavar Burton for his contribution to this code and my childhood.
# 2. Nice spacing
# 3. Move your file reading method to its own class **SEEMS UNNECESSARY**
# 

# In[ ]:




