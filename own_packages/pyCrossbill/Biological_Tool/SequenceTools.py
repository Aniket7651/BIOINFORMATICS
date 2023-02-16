#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 01:28:28 2021
here is the file for align the sequence and translate, includes:
--from text file
--needleman wunsch algorithm
--complimentary sequence
--serial comparison
--smith waterman algorithm
--translation
--message RNA
@author: ANIKET YADAV
"""
# numpy import portion 
import numpy as np


# function thats work as, to read text of a FASTA file ie. except first line it's read
# all line which is present in file
# returns joint string of the sequence 
def fromTXTfile(TextFile):
    try:
        seqString = []
        sequence = []
        stringForm = ''
        # open file with location
        with open(TextFile, 'r') as f:
            # except first line read whole sequence residue
            data = f.readlines()[1:]
            # defining sequence as a string
            seqString = data
            f.close()
            for element in seqString:
                sequence.append(element.replace('\n', '')) # replacing new line from FASTA file
                #  JOIN all the codons and return it
        return stringForm.join(sequence)
    except:
        # rether in case of an error
        return 'ERR: specify wrong path!'


# this is the method of DNA sequencing which is needleman Wunsch algorithum 
# developed by saul B. needleman and christian D. Wunsch and publish in 1970 
# -according to wikipedia
# returns numeric type value of sequence metrix only..
def NeedlemanWunschMatric(sequence1, sequence2, gap_penalty, match_score, mismatch_score):
    try:
        # nested function for match two nucleotides
        def match(a, b):
            if a == b:      # if both are same or matched so, returned score
                return match_score
            else:               # else get's mis matched score 
                return mismatch_score

        m = len(sequence1) + 1       # make two 'm'and 'n' metrix, increase size matrix by one  
        n = len(sequence2) + 1       # same for metrix another which is metrix 'n'
        # make and fill this metrix with '0' which is naming as 
        # !!!!!!!! 'row' !!!!!!!!
        row = np.zeros((m, n))   
        gap = gap_penalty    # gap for gap penalty
        for i in range(m):
            primaryScore = i * gap   # primary score multiply with gap penalty
            # initialize with first col. of first row
            row[i][0] = primaryScore
        # all things are same for matrix 'n' here,
        for i in range(n):
            primaryScore = i * gap
            row[0][i] = primaryScore
        # continue with second col. of second row
        for i in range(1, m):
            for j in range(1, n):
                # check diagonally on metrix as 'D' and add match or mismatch score
                # if there is match or not, also for match check one back of sequence on metrix 
                D = row[i - 1][j - 1] + match(sequence1[i - 1], sequence2[j - 1])
                H = row[i][j - 1] + gap   # horizontally and veritically add gap
                V = row[i - 1][j] + gap
                # find the max at last for appending on metrix
                row[i, j] = max(D, H, V)
                # returns metrix which name as 'row'
        return row
    except ValueError:
        # in case of value error
        return 'ERR: try it care fully again, with right arguments!'
    except IndexError:
        # or in case of indexing error
        return 'ERR: add large sequence as sequence1!\nsmaller sequence as sequence2!'


# finding complementary of a genetic sequence
# returns string type sequential value
def ComplementarySeq(sequence):
    try:
        complementarySeq = ''
        # data for replacing a code to another, store in dict. form
        replace = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A',
                   'a': 'T', 't': 'A', 'g': 'C', 'c': 'G'}
        for i in sequence:
            key = i         # check again and again, is 'i' in replace dict or not.. 
            if i == key:
                complementarySeq += replace[key]     # adding to empty complementarySeq string
        return complementarySeq
    except (KeyError, ValueError):
        # in the case of getting wrong letter, return err massege 
        return 'ERR: used codon should be character only (i.e. A, T, G, C)'


# for finding reverse of the complementary of the sequence
def rev_Complementary(sequence):
    try:
        complementarySeq = ''
        # data for replacing a code to another, store in dict. form
        replace = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A',
                   'a': 'T', 't': 'A', 'g': 'C', 'c': 'G'}
        for i in sequence:
            key = i         # check again and again, is 'i' in replace dict or not.. 
            if i == key:
                complementarySeq += replace[key]     # adding to empty complementarySeq string
        return complementarySeq[::-1]
    except (KeyError, ValueError):
        # in the case of getting wrong letter, return err massege 
        return 'ERR: used codon should be character only (i.e. A, T, G, C)'


class SerialComparison():
    # read sequences 
    def __init__(self, sequence1, sequence2):
        self.sequence1 = sequence1
        self.sequence2 = sequence2

    def sequence1(self):
        return self.sequence1

    def matching(self):
        match = ''
        for s1, s2 in zip(self.sequence1, self.sequence2):
            if s1 == s2:
                match += '|'
            else:
                match += ' '
        return match

    def sequence2(self):
        return self.sequence2

    def display(self):
        return f'{self.sequence1}\n{self.matching()}\n{self.sequence2}'


def SmithWaterman(sequence1, sequence2, gap_penalty, match_score, mismatch_score):
    try:
        def match(a, b):
            if a == b:
                return match_score
            else:
                return mismatch_score

        m = len(sequence1) + 1
        n = len(sequence2) + 1
        row = np.zeros((m, n))
        gap = gap_penalty
        for i in range(1, m):
            for j in range(1, n):
                D = row[i - 1][j - 1] + match(sequence1[i - 1], sequence2[j - 1])
                H = row[i][j - 1] + gap
                V = row[i - 1][j] + gap
                row[i, j] = max(D, H, V, 0)

        def traceBack(row, sequence, B_='', oldIndex=0):
            flipMatric = np.flip(np.flip(row, 0), 1)
            i_, j_ = np.unravel_index(flipMatric.argmax(), flipMatric.shape)
            i, j = np.subtract(row.shape, (i_ + 1, j_ + 1))
            if row[i, j] == 0:
                return B_, j
            B_ = sequence[j - 1] + '-' + B_ if oldIndex - i > 1 else sequence[j - 1] + B_
            return traceBack(row[0:i, 0:j], sequence, B_, i)

        return row, traceBack(row, sequence1)
    except ValueError:
        return 'ERR: try it carefully again, with specified arguments!'
    except IndexError:
        return 'ERR: add large sequence as sequence1!\nsmaller sequence as sequence2!'


def Translation(sequence, Translate_From = 0):
    try:
        aminoAcidChain = ""
        dic = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
               'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
               'TAT': 'Y', 'TAC': 'Y', 'TAA': '_', 'TAG': 'O',
               'TGT': 'C', 'TGC': 'C', 'TGA': 'U', 'TGG': 'W',

               'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
               'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
               'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

               'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
               'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
               'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
               'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',

               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
               'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
               'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
               'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

        for i in range(Translate_From, len(sequence), 3):
            amino = sequence[i:i + 3]
            if amino in dic:
                aminoAcidChain += dic[amino]
        return aminoAcidChain
    except (KeyError, ValueError, IndexError):
        return 'ERR: oops! check file path also maybe valid letters'
   

def CDS(sequence):
    cds_ = ''
    aminoAcid_count = 0
    start_codon = sequence.find('ATG')
    stop_codon = sequence[start_codon: ].find('TAA')
    for i in range(start_codon, stop_codon+2, 3):
        aminoCodon = sequence[i:i+3]
        cds_ += aminoCodon
        aminoAcid_count += 1

    return cds_, aminoAcid_count, Translation(cds_, 0)


def MessageRNA(sequence):
    mrna = ''
    for i in sequence:
        if i == 'T':
            mrna += i.replace('T', 'U')
        else:
            mrna += i
    return mrna


############################### main call #########################################

path = 'A:\Bioinformatics\ChemyCrossbill\docs\samples\dengu3_strain-D00-0107.txt'
# s1 = GGTTGACTA, s2 = TGTTACGG
s2 = 'GGTTGACTA'
#    '-\\\----'
s3 = 'TGTTACGG'
# print(NeedlemanWunschMatric(s2, s3, -2, 1, -1))
# seq = fromTXTfile(path)
# print(seq, end='\n\n')
# print(ComplementarySeq(seq))
# print(SmithWaterman(s2, s3, -2, +1, -1))
# print(traceBack(m, s3))
# print(SerialComparison(s2, s3).display())
# print(Translation(seq))
# print(MessageRNA(seq))
