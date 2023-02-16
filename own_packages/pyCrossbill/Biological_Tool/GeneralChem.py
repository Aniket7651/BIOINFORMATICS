#!/usr/bin/env python3.10
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 01:28:28 2021
general chemistry contains:
--peptide_Molecular_weight
--ssDNA_MoleculerWeight
--ssRNA_MoleculerWeight
@author: ANIKET YADAV
"""
# for find out the number of bases
# returns the four integer type output
def no_of_code(nucleotide_seq):
    # gets the parameter as DNA sequence not "aplied for RNA" 
    (no_A, no_T, no_G, no_C) = (0, 0, 0, 0)
    for i in nucleotide_seq:
        if i == 'A':
            no_A += 1
        if i == 'T':
            no_T += 1
        if i == 'G':
            no_G += 1
        if i == 'C':
            no_C += 1
    return no_A, no_T, no_G, no_C    


# returns aproxmatly weight of protein molecules 
# returns integer type output value   
def peptide_Molecular_weight(protein_seq):
    no_of_amino_acid = 0
    for i in protein_seq:
        no_of_amino_acid += 1
        # multiply with mean value of amino acid in kiloDalton 
    mw = no_of_amino_acid * 0.11
    return mw     # answer is given in "kilo Dalton"


# for finding the value of moleculer weight of single standerd of DNA
# returns the integer type of single value output   
def ssDNA_MoleculerWeight(nucleotide_seq):   # parameter for DNA
    (no_A, no_T, no_G, no_C) = no_of_code(nucleotide_seq)
    MW_ssDNA = (no_A * 313.2) + (no_T * 304.2) + (no_G * 329.2) + (no_C * 289.2) + 79.0
    return MW_ssDNA     # returns moleculer weight in "gram/mol" 


# short function for finding temperature on which melt the perticular DNA sequence 
# returns single output of integer type
def Melting_Temperature(DNA_sequence):
    # input are checked by as a DNA sequence "not for RNA"
    a, t, g, c = no_of_code(DNA_sequence)
    # eqution  returns °C answer
    return 2*a+t + 4*g+c   # returns answer in "°C"  


# finding moleculer weight of single standred RNA sequence
# returns single integer type of output of value
def ssRNA_MoleculerWeight(nucleotide_seq):
    # takes single input as mRNA sequence parameter
    (no_A, no_U, no_G, no_C) = no_of_code(nucleotide_seq)
    MW_ssRNA = (no_A * 329.2) + (no_U * 306.2) + (no_G * 345.2) + (no_C * 305.2) + 159.0
    return MW_ssRNA     # answer will apear in "gram/mole"


# here is the short code for finding the simple concentration of the double S. DNA,
# single S.DNA and single S. RNA respectivlly 
def dsDNA_concentration(dilution_factor, spectrophotometer_reading):
    # takes two parameter for every concentration finding 
    return dilution_factor * spectrophotometer_reading * 50   # returns the answer in "microGram/ML"

def ssDNA_concentration(dilution_factor, spectrophotometer_reading):
    return dilution_factor * spectrophotometer_reading * 37   # "microGram/ML"

def ssRNA_concentration(dilution_factor, spectrophotometer_reading):
    return dilution_factor * spectrophotometer_reading * 40   # "microGram/ML"


# the function for finding the avg. (mean) of base pairs in DNA
# returns the four integer type value of output
def mean_of_basePairs(sequence):
    # input as a single DNA sequence exists
    Total_ = 0
    (countA, countT, countG, countC) = no_of_code(sequence)
    for i in sequence:
        Total_ += 1
    # returns the avg. of adinine, thymine, guanine, cytosine respectivly four output
    return countA/Total_, countT/Total_, countG/Total_, countC/Total_


# for read emprical formula of a chemical structure, only compatible on organic compound
# to get input as a emprical formula of drug and gives output as a two list:-
# first is number of element which is contain in formula, and second one is name of the element  
def read_empirical(empirical_formula):
    import re
    # re is the regular expression for finding numeric in emprical
    # empirical like this = 'C10H16N2O2'
    element = []
    # finding the organic element into emprical 
    for i in empirical_formula:
        if i == 'C':            # find carbon
            element.append(i)
        elif i == 'H':          # find hydorgen
            element.append(i)
        elif i == 'N':          # find nitrogen
            element.append(i)
        elif i == 'O':          # find oxygen
            element.append(i)
        elif i == 'F':          # find florine
            element.append(i)
        elif i == 'I':          # find iodine
            element.append(i)
    # sepratly find bromine, clorine  
    if 'Br' in empirical_formula:
        element.append('Br')
    if 'Cl' in empirical_formula:
        element.append('Cl')
    # find how many element present in a particular compound 
    no_contant_list = re.findall(r"\d+", empirical_formula)
    # return number of element and name of element
    return no_contant_list, element


# block of code which find molecular weight using emprical formula('read_empirical' function above)
# it's take input as return function of read_empirical and you gets output a integer type mol. weight 
def Organic_MolecularWeight(no_contant_list, element):
    no_cont = 0
    no_element = 0
    # count number of containt and element 
    # for genrating a condition, is both are equal -> so find the mol. weight at the line 134.
    for a in no_contant_list:
        no_cont += 1
    for b in element:
        no_element += 1
    if no_cont == no_element:
        periodic = {'H': 1.007, 'C': 12.01, 'N': 14.006, 'O': 15.999, 'S': 32.065,
                    'F': 18.998, 'Cl': 35.453, 'Br': 79.904, 'I': 126.90}
        # there is the mol. weight of a perticular organic element
        Mw = 0
        # initially mol. weight is zero
        for ele, c in zip(element, no_contant_list):
            if ele in periodic:
                # multiply element, if it's present in periodic dictionary 
                # and perticular number of element
                Mw += periodic[ele]*int(c)
        # return mol. weight
        return Mw
    # if both are not equal -> so add 1 on number of element list type
    # means only single element is present, not multiple form
    else:
        diff = no_element - no_cont
        for n in range(0, diff):
            no_contant_list.append(1)
        # at the last recall all over function(by recursion)
        return Organic_MolecularWeight(no_contant_list, element)


# this is the program for predict the value of boiling point chemical compound using joback method
# group contribution is taken from wikipidia site, value of group contribution
# is acually the value of phase transtions of tempratures, use only on organic compound 
def boiling_group_contribution(smile_represnt):
    # takes the smile represent of chemical structure
    # SMALL letter for ring group and CAPITAL for non-ring 

    forCarbon = {'C': 22.96, 'c': 22.96}
    forSulfur = {'S': 61.48, 's': 61.48}
    forHalogen = {'F': -0.03, 'I': 93.84}
    forNitrogen = {'N': 66.10, 'n': 64.48}
    forOxygen = {'O': 70.65, 'o': 70.65}
    # this all value is the mean of Tb (phase transtions) for acordinglly ringed or non-ringed group
    sum_Gi = 0
    if ('Cl'or'Br' or 'NOO') in smile_represnt:
        dobCounter = (smile_represnt.count('Br')*66.86,
                    smile_represnt.count('Cl')*38.13,
                    smile_represnt.count('NOO')*152.54)
        # for Cl, Br, NOO is use sepretly
        # also add in sum of Gi (group contribution)
        for num in dobCounter:
            sum_Gi += num
    for group in smile_represnt:
        if group in forCarbon:
            sum_Gi += forCarbon[group]   # add corbon Tb in sum of Gi
        elif group in forHalogen:       
            sum_Gi += forHalogen[group]     # add halogen Tb in sum of Gi
        elif group in forNitrogen:
            sum_Gi += forNitrogen[group]        # add nitrogen Tb in sum of Gi
        elif group in forSulfur:
            sum_Gi += forSulfur[group]              # add sulfur Tb in sum of Gi
        elif group in forOxygen:
            sum_Gi += forOxygen[group]          # add oxygen Tb in sum of Gi
    # returns answer in floating point use in prediction of normal boiling point 
    return sum_Gi

# NBP is stands for normal boiling point using joback method
# returns single integer type of value 
def NBP(smile_represnt):
    # takes smile represent of chemical structure
    Tb_K = 198.20 + boiling_group_contribution(smile_represnt)
    # returns the output into kelvin
    return Tb_K


# this is the program for predict the value of melting point chemical compound using joback method
# group contribution is taken from wikipidia site, value of group contribution
# is acually the value of phase transtions of tempratures, use only on organic compound 
def melting_group_contribution(smile_represnt):
    # takes the smile represent of chemical structure
    # SMALL letter for ring group and CAPITAL for non-ring 

    forCarbon = {'C': 18.97, 'c': 18.97}
    forSulfur = {'S': 44.80, 's': 44.80}
    forHalogen = {'F': -15.78, 'I': 41.69}
    forNitrogen = {'N': 66.67, 'n': 66.67}
    forOxygen = {'O': 55.78, 'o': 55.78}
    # this all value is the mean of Tb (phase transtions) for acordinglly ringed or non-ringed group
    sum_Gi = 0
    if ('Cl'or'Br' or 'NOO') in smile_represnt:
        dobCounter = (smile_represnt.count('Br')*43.43,
                    smile_represnt.count('Cl')*13.55,
                    smile_represnt.count('NOO')*127.24)
        for num in dobCounter:
            sum_Gi += num
        # for Cl, Br, NOO is use sepretly
        # also add in sum of Gi (group contribution)
        # # for Cl, Br, NOO is use sepretly
    for group in smile_represnt:
        if group in forCarbon:
            sum_Gi += forCarbon[group]      # add corbon Tb in sum of Gi
        elif group in forHalogen:
            sum_Gi += forHalogen[group]         # add halogen Tb in sum of Gi
        elif group in forNitrogen:
            sum_Gi += forNitrogen[group]        # add nitrogen Tb in sum of Gi
        elif group in forSulfur:
            sum_Gi += forSulfur[group]          # add sulfur Tb in sum of Gi
        elif group in forOxygen:
            sum_Gi += forOxygen[group]          # add oxygen Tb in sum of Gi
    # returns answer in floating point use in prediction of normal melting point
    return sum_Gi

# NBP is stands for normal melting point using joback method
# returns single integer type of value 
def NMP(smile_represnt):
    Tm_K = 122.5 + melting_group_contribution(smile_represnt)
    # returns the output into kelvin
    return Tm_K


# problem is there conversion in temperature
# using simple formula of temp. conversion
class Temperature_Conversion():
    def celsius_to_kelvin(self, celsius_):
        # return answer in K
        return celsius_ + 273.15
    
    def kelvin_to_celsius(self, kelvin_):
        # return answer in °C
        return kelvin_ - 273.15

    def fahrenheit_to_celsius(self, fahrenheit_):
        # # return answer in °C
        return fahrenheit_ - 32 * 0.55
    
    def celsius_to_fahrenheit(self, celsius_):
        # return answer in °F
        return celsius_ * 0.55 + 32

    def fahrenheit_to_kelvin(self, fahrenheit_):
        # # return answer in K
        return fahrenheit_ - 32 * 0.55 + 273.15
    
    def kelvin_to_fahrenheit(self, kelvin_):
        # return answer in °F
        return kelvin_ - 273.15 * 0.55 + 32


########################### END OF THE PROGRAM ################################

# print(peptide_Molecular_weight(s2))
# print(ssDNA_MoleculerWeight(s2))      # 3465.2 g/mol
# print(ssRNA_MoleculerWeight(s1))
# no_contant_list, element = read_empirical('C16H25NO2')
# print('MW: ', Organic_MolecularWeight(no_contant_list, element))
smile_represnt = 'CCOC(=O)C1=CC2CN(C(=O)C3=C(N2C=N1)C=CC(=C3)F)C' #  ASIPRINE
# CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3  DIAZEPAM
# 2-[2-(2,6-dichloroanilino)phenyl]acetic acid -- diclofenac
# print('boiling: ',NBP(smile_represnt))
# print('melting: ', NMP(smile_represnt))
# print(dsDNA_concentration(50, 0.65))
# print(ssRNA_concentration(50, 0.65))
# print(ssDNA_concentration(50, 0.65))
# print(Melting_Temperature(s2))
# print(mean_of_basePairs(s2))
# print(Temperature_Conversion().kelvin_to_celsius(458.14))