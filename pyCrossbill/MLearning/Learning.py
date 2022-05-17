#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 01:28:28 2021

@author: ANIKET YADAV
"""


def calculate_error(pred_list, actual_list):
    errors = []
    for pred_, actual in zip(pred_list, actual_list):
        errors.append(actual-pred_)
    return errors

# return two values which is more then 100 and less then 100
def validate_error(errors):
    (less_then100, more_then100) = ([], [])
    for i in errors:
        if i < 100:
            less_then100.append(i)
        else:
            more_then100.append(i)
    return less_then100, more_then100


################################## END OF THE PROGRAM ###############################################
lists = [100, 29, 173, 728, 26, 93, 103]
print(validate_error(lists))