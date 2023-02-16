#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 19:08:16 2021
create your own file formats, includes:
--from txt file
--fasta file
@author: ANIKET YADAV
"""

import os

# this program for read the fasta file in python 
# returns string type of single lined sequence
def readFASTAfile(TextFile):
    try:
        seqString = []
        sequence = []
        stringForm = ''
        with open(TextFile, 'r') as f:
            # read from second line, drop first description line
            data = f.readlines()[1:]
            seqString = data
            f.close()
            for element in seqString:
                # repplace new line to serice
                sequence.append(element.replace('\n', ''))
        # returns a string type in a single line
        return stringForm.join(sequence)
    # if any type of error in reading 
    except:
        return 'ERR: specify wrong path!'

# create fasta file using locus, DEFINITION, ORGANISM, VERSION and file_Path of text file of
# sequence return a new file where you give file path
def makeFASTAFile(sequence, LOCUS, DEFINITION, ORGANISM, VERSION, file_Path):
    try:
        with open(file_Path, 'a') as f1:
            # write first description line
            f1.write(f'>{LOCUS}| {VERSION}| {ORGANISM}| {DEFINITION}\n')
            count = 0
            # read a item and add into a single line
            for i in sequence:
                f1.write(i)
                count += 1
                # check the condition if it's item count (codon) is 65
                # then add a new line also reset the count number to 0 
                if count == 65:
                    f1.write('\n')
                    count = 0
        # returns the path of the given file where save the created fasta file 
        return file_Path
    except:
        return "somthing went's wrong here, please try again!"



################################ END OF THE PROGRAM ########################################

path = os.path.abspath('pyCrossbill/docs/samples/dengu3_strain-D00-0107.txt')
# sequence = fromTXTfile(path)
print(readFASTAfile(path))
# print(FASTAFile(sequence, 'NC_500.64', 'dengu3_strain-D00-0107', 'dengu3', '546.6', 'file.fasta'))
