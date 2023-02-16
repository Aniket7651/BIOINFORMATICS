#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 01:30:21 2021
this is measure the distance between two sequence, includes:
--lempele ziv
--lempel ziv complexity
--euclidean distance
@author: ANIKET YADAV
"""
import numpy as np


def LempeleZiv(sequence):
    key_dic = {}
    index = 0
    increment = 1
    countElement = 0

    while True:
        if not (len(sequence) >= index + increment):
            break
        subString = sequence[index: index + increment]
        # print(subString, index, increment)
        if subString in key_dic:
            increment += 1

        else:
            key_dic[subString] = 0
            index += increment
            increment = 1
    listOflemple = list(key_dic)
    for element in listOflemple:
        if element:
            countElement += 1

    return listOflemple, countElement


def NormalizeCompressionDistance(LempeleZiv_of_seq1, LempeleZiv_of_seq2, both_LempleZiv):
    C_x = LempeleZiv_of_seq1
    C_y = LempeleZiv_of_seq2
    C_xy = both_LempleZiv

    numrator = C_xy - min(C_x, C_y)
    NCD = numrator / max(C_x, C_y)

    return NCD


def EuclideanDistance(sequence1, sequence2):
    def words(sequence):
        lis1, lis2 = [], []
        for i in range(0, len(sequence), 3):
            lis1.append(sequence[i:i + 3])
        for i in range(1, len(sequence), 3):
            lis2.append(sequence[i:i + 3])

        def counter(lis):
            count = 0
            for i in lis[-1]:
                count += 1
            if count == 3:
                validList = lis[:]
            else:
                validList = lis[:-1]
            return validList

        nlis1, nlis2 = counter(lis1), counter(lis2)
        return nlis1 + nlis2

    W_X = words(sequence1)
    W_Y = words(sequence2)
    W_x, W_y = set(W_X), set(W_Y)
    union_of_set = W_x.union(W_y)
    cx, cy = [], []
    for i in union_of_set:
        cx.append(W_X.count(i))
        cy.append(W_Y.count(i))
    diff = []
    for i, j in zip(cx, cy):
        d = i - j
        diff.append(d ** 2)
    euclidean = (np.sum(diff)) ** 0.5
    print(W_X, W_Y)
    print('union: ', union_of_set)
    return euclidean


################################# main call ####################################

seq1 = 'CATGTG'
seq2 = 'CATGTT'
both = 'CATGTGCATGTT'
print(LempeleZiv(both))
# print(NormalizeCompressionDistance(5, 5, 7))
print(EuclideanDistance(seq1, seq2))
