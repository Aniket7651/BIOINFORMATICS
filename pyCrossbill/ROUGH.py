from GeneralChem import (ssRNA_MoleculerWeight, mean_of_basePairs)

def fromTXTfile(TextFile):
    try:
        seqString = []
        sequence = []
        stringForm = ''
        with open(TextFile, 'r') as f:
            data = f.readlines()[1:]
            seqString = data
            f.close()
            for element in seqString:
                sequence.append(element.replace('\n', ''))
        return stringForm.join(sequence)
    except:
        return 'ERR: specify wrong path!'


class sequenceInfo():
    def __init__(self, fastaFile):
        self.fastaFile = fastaFile 

    def display(self):
        txtseq = fromTXTfile(self.fastaFile)
        return f'{ssRNA_MoleculerWeight(txtseq)}\n{mean_of_basePairs(txtseq)}'

        
def melting_group_contribution(smile_represnt):
    # SMALL letter for ring group and CAPITAL for non-ring 

    forCarbon = {'C': 15.17, 'c': 26.58}
    forSulfur = {'S': 27.24, 's': 79.93}
    forHalogen = {'F': -15.78, 'I': 41.69}
    forNitrogen = {'N': 59.43, 'n': 84.95}
    forOxygen = {'O': 62.31, 'o': 60.61}
    sum_Gi = 0
    if ('Cl'or'Br' or 'NOO') in smile_represnt:
        dobCounter = (smile_represnt.count('Br')*43.43,
                    smile_represnt.count('Cl')*13.55,
                    smile_represnt.count('NOO')*127.24)
        for num in dobCounter:
            sum_Gi += num
    for group in smile_represnt:
        if group in forCarbon:
            sum_Gi += forCarbon[group]
        elif group in forHalogen:
            sum_Gi += forHalogen[group]
        elif group in forNitrogen:
            sum_Gi += forNitrogen[group]
        elif group in forSulfur:
            sum_Gi += forSulfur[group]
        elif group in forOxygen:
            sum_Gi += forOxygen[group]
    return sum_Gi/2


def NMP(smile_represnt):
    Tm_K = 122.5 + melting_group_contribution(smile_represnt)
    # returns the output into kelvin
    return Tm_K


class Temperature_Conversion():
    def celsius_to_kelvin(self, celsius_):
        return f'{celsius_ + 273.15} K'
    
    def kelvin_to_celsius(self, kelvin_):
        return f'{kelvin_ - 273.15} °C'

    def fahrenheit_to_celsius(self, fahrenheit_):
        return f'{fahrenheit_ - 32 * 0.55} °C'
    
    def celsius_to_fahrenheit(self, celsius_):
        return f'{celsius_ * 0.55 + 32} °F'

    def fahrenheit_to_kelvin(self, fahrenheit_):
        return f'{fahrenheit_ - 32 * 0.55 + 273.15} K'
    
    def kelvin_to_fahrenheit(self, kelvin_):
        return f'{kelvin_ - 273.15 * 0.55 + 32} °F'


def read_empirical(empirical_formula):
    import re
    # empirical like this = 'C10H16N2O2'
    element = []
    for i in empirical_formula:
        if i == 'C':
            element.append(i)
        elif i == 'H':
            element.append(i)
        elif i == 'N':
            element.append(i)
        elif i == 'O':
            element.append(i)
        elif i == 'F':
            element.append(i)
        elif i == 'I':
            element.append(i)
    if 'Br' in empirical_formula:
        element.append('Br')
    if 'Cl' in empirical_formula:
        element.append('Cl')
    no_contant_list = re.findall(r"\d+", empirical_formula)
    return no_contant_list, element

def Identify_IUPAC(IUPAC):
    # 3-cloro-3-methylhexane
    bond_tokens = ['ane', 'ene', 'yne']
    carbon_tokens = ['meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non', 'dec']

from PIL import Image
from numpy import asarray
image = Image.open('A:/Bioinformatics/ChemyCrossbill/docs/tramadol.png')
arr_img = asarray(image)
print(arr_img.shape)
path2 = 'A:/Bioinformatics/ChemyCrossbill/test.txt'
# seq = fromTXTfile(path)
smile = 'c1ccccn1'
# no_contant_list, element = read_empirical('C10H16N2O2')
# print(Organic_MolecularWeight(no_contant_list, element)) 
# print(NMP(smile), 'K')
# print(Temperature_Conversion().kelvin_to_celsius(NMP(smile)))
# print(sequenceInfo(path2).display())