#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 01:28:28 2021
contains machine learning programs for predicting the biological data using algorithms 
@author: ANIKET YADAV
"""

# nessesary imports
import csv, os

# pre loaded dataset path which use as default dataset file 
# dataset = os.path.abspath('pyCrossbill/MLearning/drugset.csv')


# block of program for reading csv file
class load_csv():
    """
    For reading a simple CSV file by calling load_csv(file path).dataset(select_column, select_row) method\n
    --load_csv(file path).top(print_upto_row)\n
    --load_csv(file path).tail(print_upto_row)\n
    """
    # specifying path by init
    def __init__(self, file_path):
        self.file = file_path

    # display dataset by two arguments column and upto
    # column = all ----- for represent all the column take's only all, and other's integer
    # upto = all ------- for representing all the data of a csv file, row also take's only all, and other's integer
    def dataset(self, column='all', upto='all'):
        '''Function of class load_csv() to use for represent the CSV data\n
        column='all'-- is select all column by default, you can select number, how many column you wish to select\n
        upto='all'-- is select all row by default, you can select number, like how many rows you wish to select'''
        try:
            # step to store data in row list variable
            rows = []
            # open file and manipulate as csvf instance
            with open(self.file, 'r') as csvf:
                for row in csv.reader(csvf):
                    # check if all column is set 
                    # print all row 
                    if column == 'all':
                        rows.append(row[:])
                    # else print a single particular column 
                    else:
                        rows.append(row[column])
                # close the csv file
                csvf.close()
            # check if upto selected all
            # it's return all row
            if upto == 'all':
                return rows[1:]
            else:
                # else return row upto, where you are selected
                return rows[1:upto+1]
        except IndexError:
            # return error if you are selected wrong column
            return f"column {column} not exist in your dataset"
     
    def top(self, upto):
        '''top(upto_row) method to select row upto you wish to select from 0th row'''
        # firstly select all row and column
        select = self.dataset()
        # select dataset from zero to row upto....
        return select[:upto]

    def tail(self, upto):
        '''tail(upto_row) method to select row upto you wish to select from last data row'''
        select = self.dataset()
        # same as for tail select vary last occurence in data row and move last to upward side....
        return select[-upto:]


class loss_function():
    """
    Loss function is determine the accuracy of our
    """
    def __init__(self, predicted_list, actual_list):
        self.predicted = predicted_list
        self.actual = actual_list
    
    # mean square error (MSE)
    def L2_loss(self):
        error = 0.0
        N = len(self.actual)
        for pred_, actual in zip(self.predicted, self.actual):
            error += (pred_ - actual)**2
        return error/float(N)

    def absolute_error(self):
        errors = []
        for pred_, actual in zip(self.predicted, self.actual):
            errors.append(pred_ - actual)
        return errors
    
    # mean absolute error (MAE)
    def L1_loss(self):
        errors = 0.0
        N = len(self.actual)
        for pred_, actual in zip(self.predicted, self.actual):
            errors += pred_ - actual
        return errors/float(N)
   
    def SSE(self):
        errors = []
        for pred_, actual in zip(self.predicted, self.actual):
            errors.append((pred_ - actual)**2)
        return errors

    def RMSD(self):
        N = len(self.actual)
        return (self.L2_loss()/N)**0.5

    def hinge(self):
        hinge = []
        for t, y in zip(self.actual, self.predicted):
            hinge.append(max(0, 1 - t * y))
        return hinge


class statistics():
    def __init__(self, x_set):
        self.x = x_set
        self.N = len(x_set)

    def mean(self):
        Elements = sum(self.x)
        return Elements/float(self.N)

    def variance(self):
        mean = self.mean()
        sumation = 0.0
        for i in self.x:
            sumation += (i - mean)**2
        return sumation/self.N

    def median(self):
        if self.N % 2 == 0:
            term = int((self.N + 1)/2)
            return self.x[term-1]
        else:
            d = self.N/2
            term = int(d + (d +1)/2)
            return self.x[term-1]    
    
    def standard_dev(self):
        mean = self.mean()
        sumation = 0.0
        for i in self.x:
            sumation += (i - mean)**2
        return (sumation/self.N)**0.5


def covarince(x_set, y_set):
    x_mean = statistics(x_set).mean()
    y_mean = statistics(y_set).mean()
    addtion = 0.0
    n = len(x_set)
    for x, y in zip(x_set, y_set):
        addtion += (x - x_mean)*(y - y_mean)
    return addtion/n


# we need a way to determine if thhere is linear correlation or not, so we calculate what is know 
# as the PRODUCT-MOMENT CORRELATION COEFFICIENT.
def Product_moment_CC(x_set, y_set):
    std_y = statistics(y_set).standard_dev()
    std_x = statistics(x_set).standard_dev()
    x_mean = statistics(x_set).mean()
    y_mean = statistics(y_set).mean()
    addtion = 0.0
    n = len(x_set)
    for x, y in zip(x_set, y_set):
        addtion += (x - x_mean)*(y - y_mean)
    covar = addtion/n
    return covar/std_x*std_y


################################## END OF THE PROGRAM ###############################################
file = 'A:/BIOINFORMAICS/own_packages/pyCrossbill/datasets/drugset.csv'
pd = [12.3, 53.2, 52.2, 13.4, 83.5]
at = [26.7, 34.2, 52.0, 63.5, 12.6]
a = [0,1,1,0,1,1,1,0,0,1]
p = [0,0,0,1,0,1,1,0,1,0]
set = [3, 6, 9, 2, 7]
# print(loss_function(p, a).hinge())
# print(covarince(at, pd))
# print(load_csv(file).dataset(column=1, upto=5))
print(covarince(at, pd))
# print(statistics(set).mean())