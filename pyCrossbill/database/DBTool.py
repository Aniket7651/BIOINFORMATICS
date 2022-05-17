#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 01:41:36 2021

@author: linux_kali
"""

import sqlite3


def InsertJob(locus, definition, version, organism, refrence, authors, comments, origin, database='NucleotideDB.db'):
    try:
        conn = sqlite3.connect(database)
        cur = conn.cursor()
        query = """INSERT INTO Nucleotide (LOCUS, DEFINITION, VERSION, ORGANISM,
                Ref, AUTHORS, COMMENT, ORIGIN) VALUES (?,?,?,?,?,?,?,?);"""
        data = (locus, definition, version, organism, refrence, authors, comments, origin)
        cur.execute(query, data)
        conn.commit()
        cur.close()
        return f'successfully {database}'
    except sqlite3.Error as err:
        return f'insertion Error {err}'
    # finally:
    #    if conn:
    #        conn.close()
    #        return 'closed'


database = '/home/linux_kali/Documents/Bioinformatics/PyBio-Information/DataBases/NucleotideDB.db'
print(InsertJob('NC_536.45', 'this is virus', 54, 'virus', 7.5, 'ani', 'affg', 'gcatgcagctagctagctacgtagcatc'))
