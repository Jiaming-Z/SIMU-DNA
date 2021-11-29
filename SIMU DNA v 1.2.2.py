#!/usr/bin/env python
# coding: utf-8

# In[20]:


# phase 1: protein synthesis!
DNA_5to3 = ["A", "T", "C", "G", "A", "A", "A", "C", "C", "T", "G", "G"]
def protein_synthesis_main(DNA_5to3):
    full_double_helix(DNA_5to3)
    pre_mRNA = transcription(DNA_5to3)
    mRNA = RNA_splicing(pre_mRNA, intron_indices)
    mRNA_codons = make_codon(mRNA)
    tRNA = find_tRNA(mRNA_codons)
    return peptide_forming(mRNA_codons)
protein_synthesis_main(DNA_5to3)


# In[5]:


DNA_5to3 = ["A", "T", "C", "G", "A", "A", "A", "C", "C", "T", "G", "G"]
def full_double_helix(DNA_5to3):
    DNA_3to5 = []
    full_DNA = []
    for i in range(len(DNA_5to3)):
        if DNA_5to3[i] == "A":
            DNA_3to5.append("T")
        elif DNA_5to3[i] == "T":
            DNA_3to5.append("A")
        elif DNA_5to3[i] == "C":
            DNA_3to5.append("G")
        elif DNA_5to3[i] == "G":
            DNA_3to5.append("C")
        full_DNA.append([DNA_5to3[i], DNA_3to5[i]])
    print("lagging strand:", DNA_3to5)
    print("full double helix:", full_DNA)
full_double_helix(DNA_5to3)


# In[6]:


DNA_5to3 = ["A", "T", "C", "G", "A", "A", "A", "C", "C", "T", "G", "G"]
def transcription(DNA_5to3):
    pre_mRNA = []
    for i in range(len(DNA_5to3)):
        if DNA_5to3[i] == "A":
            pre_mRNA.append("U")
        elif DNA_5to3[i] == "T":
            pre_mRNA.append("A")
        elif DNA_5to3[i] == "C":
            pre_mRNA.append("G")
        elif DNA_5to3[i] == "G":
            pre_mRNA.append("C")
    return pre_mRNA
print("pre_mRNA:", transcription(DNA_5to3))


# In[7]:


pre_mRNA = transcription(DNA_5to3)
intron_indices = [[0, 3], [4, 6], [9, 10]] #deletion will include first index, exclude second index

def RNA_splicing(pre_mRNA, intron_indices):
    new_copy_mRNA = pre_mRNA[:]
    frameshift_lag = 0
    for i in range(len(intron_indices)):        
        if i == 0:
            new_copy_mRNA = new_copy_mRNA[:intron_indices[i][0]] + new_copy_mRNA[intron_indices[i][1]:]
        else:
            frameshift_lag += intron_indices[i-1][1] - intron_indices[i-1][0]
            new_copy_mRNA = new_copy_mRNA[:intron_indices[i][0]-frameshift_lag] + new_copy_mRNA[intron_indices[i][1]-frameshift_lag:]
    return new_copy_mRNA
print("mRNA:", RNA_splicing(pre_mRNA, intron_indices))
    


# In[8]:


mRNA = ['C', 'U', 'G', 'G', 'C', 'C', 'G', 'C']
def make_codon(mRNA):
    mRNA_codons = []
    for i in range(0, len(mRNA), 3):
        if i+2 < len(mRNA):
            mRNA_codons.append([mRNA[i], mRNA[i+1], mRNA[i+2]])
    return mRNA_codons
make_codon(mRNA)


# In[9]:


mRNA_codons = [['C', 'U', 'G'], ['G', 'C', 'C']]
def find_tRNA(mRNA_codons):
    tRNA_chain = []    
    for codon in mRNA_codons:
        tRNA = []
        for i in range(len(codon)):
            if codon[i] == "A":
                tRNA.append("U")
            elif codon[i] == "U":
                tRNA.append("A")
            elif codon[i] == "C":
                tRNA.append("G")
            elif codon[i] == "G":
                tRNA.append("C")
        tRNA_chain.append(tRNA)
    return tRNA_chain
find_tRNA(mRNA_codons)
        


# In[10]:


def peptide_forming(mRNA_codons):
    peptide_chain = []
    for codon in mRNA_codons:
        if codon == ['U', 'U', 'U'] or codon == ['U', 'U', 'C']:
            peptide_chain.append('phenylalanine')
        elif codon == ['U', 'U', 'A'] or codon == ['U', 'U', 'G'] or codon[:2] == ['C', 'U']:
            peptide_chain.append('leucine')
        elif codon == ['A', 'G', 'C'] or codon == ['A', 'G', 'U'] or codon[:2] == ['U', 'C']:
            peptide_chain.append('serine')
        elif codon[:2] == ['C', 'C']:
            peptide_chain.append('proline')
        elif codon[:2] == ['A', 'C']:
            peptide_chain.append('threonine')
        elif codon[:2] == ['G', 'C']:
            peptide_chain.append('alanine')
        elif codon[:2] == ['G', 'U']:
            peptide_chain.append('valine')
        elif codon == ['A', 'G', 'A'] or codon == ['A', 'G', 'G'] or codon[:2] == ['C', 'G']:
            peptide_chain.append('arginine')
        elif codon[:2] == ['G', 'G']:
            peptide_chain.append('glycine')
        elif codon == ['U', 'G', 'C'] or codon == ['U', 'G', 'U']:
            peptide_chain.append('cysteine')
        elif codon == ['U', 'A', 'U'] or codon == ['U', 'A', 'C']:
            peptide_chain.append('tyrosine')
        elif codon == ['A', 'A', 'U'] or codon == ['A', 'A', 'C']:
            peptide_chain.append('asparagine')
        elif codon == ['A', 'A', 'A'] or codon == ['A', 'A', 'G']:
            peptide_chain.append('lysine')
        elif codon == ['C', 'A', 'A'] or codon == ['C', 'A', 'G']:
            peptide_chain.append('glutamine')
        elif codon == ['C', 'A', 'U'] or codon == ['C', 'A', 'C']:
            peptide_chain.append('histidine')
        elif codon == ['A', 'U', 'U'] or codon == ['A', 'U', 'C'] or codon == ['A', 'U', 'A']:
            peptide_chain.append('isoleucine')
        elif codon == ['G', 'A', 'U'] or codon == ['G', 'A', 'C']:
            peptide_chain.append('aspartic acid')
        elif codon == ['G', 'A', 'G'] or codon == ['G', 'A', 'A']:
            peptide_chain.append('glutamic acid')
        elif codon == ['U', 'G', 'G']:
            peptide_chain.append('tryptophan')
        elif codon == ['A', 'U', 'G']:
            peptide_chain.append('methionine')
        elif codon == ['U', 'A', 'A'] or codon == ['U', 'A', 'G'] or codon == ['U', 'G', 'A']:
            peptide_chain.append('STOP')
            return peptide_chain
    return peptide_chain
        
peptide_forming([['G', 'A', 'C'], ['U', 'A', 'A'], ['C', 'G', 'G']])            


# In[11]:


def count_frequencies(DNA_5to3):
    numA, numT, numC, numG = 0, 0, 0, 0
    for base in DNA_5to3:
        if base == 'A':
            numA += 1
        if base == 'T':
            numT += 1
        if base == 'C':
            numC += 1
        if base == 'G':
            numG += 1
    return [['A', 'T', 'C', 'G'], [numA, numT, numC, numG]]
count_frequencies(DNA_5to3)


# In[12]:


# phase 2: muations!


# In[13]:


def point_mutation(DNA_5to3, pos, new_base):
    new_DNA_5to3 = DNA_5to3[:]
    new_DNA_5to3[pos] = new_base
    return new_DNA_5to3
point_mutation(['A', 'T', 'C'], 1, 'A')


# In[14]:


def frameshift_insertion(DNA_5to3, pos, new_strand):
    new_DNA_5to3 = DNA_5to3[:pos] + new_strand + DNA_5to3[pos:]
    return new_DNA_5to3
frameshift_insertion(['A', 'T', 'C'], 1, ['A', 'G'])


# In[15]:


def frameshift_deletion(DNA_5to3, pos, length):
    assert length <= len(DNA_5to3[pos:])
    new_DNA_5to3 = DNA_5to3[:pos] + DNA_5to3[pos + length:]
    return new_DNA_5to3
frameshift_deletion(['A', 'T', 'C'], 1, 1)


# In[16]:


import random
def random_mutation(DNA_5to3):
    mutation_type = random.randint(0, 1)
    pos = random.randint(0, len(DNA_5to3)-1)
    if mutation_type == 0: #point
        new_base = random.choice(['A', 'T', 'C', 'G'])
        return point_mutation(DNA_5to3, pos, new_base)
    else: #frameshift
        insert_or_delete = random.randint(0, 1)
        if insert_or_delete == 0: #insertion
            insert_length = random.randint(0, 51)#upper limit 51 bases
            new_strand = [random.choice(['A', 'T', 'C', 'G']) for i in range(insert_length)]
            return frameshift_insertion(DNA_5to3, pos, new_strand)
        else: #deletion
            delete_length = random.randint(0, len(DNA_5to3)-1)
            return frameshift_deletion(DNA_5to3, pos, delete_length)
    
random_mutation(['A', 'T', 'C'])


# In[17]:


# phase 3: minimum mutations!


# In[18]:


#Wagner-Fischer Algorithm of Lenvenshtein distance, find the minimum number of single base 
#muations(insertion, deletion, substitution) required to reach from reference to target
import numpy as np 
# def minimum_mutations(reference, target):
#     dis = np.zeros((len(reference)+1, len(target)+1)) #array
#     for i in range(0, len(reference)+1):
#         dis[i, 0] = i 
#     for j in range(0, len(target)+1):
#         dis[0, j] = j
#     for j in range(1, len(target)+1):
#         for i in range(1, len(reference)+1):
#             if target[j-1] == reference[i-1]:
#                 substitution_cost = 0
#             else: 
#                 substitution_cost = 1
#             dis[i, j] = min(dis[i-1,j] + 1, dis[i,j-1] + 1, dis[i-1,j-1] + substitution_cost) # deletion, insertion, substitution respectively 
#     return dis[len(reference), len(target)]

# two-line array formating
#https://en.wikipedia.org/wiki/Levenshtein_distance#Iterative_with_two_matrix_rows
def minimum_mutations(reference, target):
    v0 = np.zeros(len(target)+2) #vector/list for previous row
    v1 = np.zeros(len(target)+2) #vector/list for current row
    for i in range(0, len(target)+2):
        v0[i]= i #previous row distance,start with empty reference
    for i in range(0, len(reference)):
        v1[0] =i+1
        for j in range(0, len(target)):
            deletion_cost =v0[j+1] + 1
            insertion_cost = v1[j] + 1
            if target[j] == reference[i]:
                substitution_cost = v0[j]
            else: 
                substitution_cost = v0[j] + 1
            v1[j+1]= min(deletion_cost, insertion_cost, substitution_cost)
        v0, v1 = v1, v0
    return v0[len(target)]


# In[19]:


minimum_mutations(DNA_5to3, ['A', 'G', 'C', 'T'])


# In[ ]:


#widgets!!


# In[ ]:





# In[ ]:




