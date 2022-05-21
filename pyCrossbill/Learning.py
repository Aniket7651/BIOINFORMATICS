#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 01:28:28 2021
contains machine learning programs for predicting the biological data using algorithms 
@author: ANIKET YADAV
"""

from GeneralChem import boiling_group_contribution, NBP, melting_group_contribution, NMP
import csv, os


def load_csv(file_path):
    rows = []
    with open(file_path, 'r') as csvf:
        for row in csv.reader(csvf):
            rows.append(row)
        csvf.close()
    return rows[1:]


def Mean_Square_error(pred_list, actual_list):
    error = 0
    N = len(actual_list)
    for pred_, actual in zip(pred_list, actual_list):
        error += (actual-pred_)**2
    return error/N


# return two values which is more then 100 and less then 100
def validate_error(errors):
    (less_then100, more_then100) = ([], [])
    for i in errors:
        if i < 100:
            less_then100.append(i)
        else:
            more_then100.append(i)
    return less_then100, more_then100


def predict_sets(dataset=os.path.abspath('pyCrossbill/MLearning/raw_drug.csv')):
    high_err = 0
    low_err = 0
    Tm_ = []
    for data in load_csv(dataset):
        g_i = melting_group_contribution(data, high_err, low_err)
        Tm_.append(122.5 + melting_group_contribution(g_i))
    



################################## END OF THE PROGRAM ###############################################
file = 'A:/BIOINFORMAICS/pyCrossbill/MLearning/raw_drug.csv'
print(load_csv(file))
# lists = [100, 29, 173, 728, 26, 93, 103]
# print(validate_error(lists))