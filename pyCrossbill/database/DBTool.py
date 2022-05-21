#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 01:41:36 2021
database programs for storing fetching and analyse protein and nucleotide data
what's you are store
@author: linux_kali
"""

# By default pkg. use sqlite3 database 
import sqlite3
import os         # for defining path of the database file 

# specifying database ________
database = os.path.abspath('pyCrossbill/database/pyCrossbill.sqlite')

# insert data to nucleotide database table
# returns about the succses massege and path of the database file   
def InsertJob_nucleotide(locus, definition, version, organism, refrence, authors, comments, origin, 
    database=os.path.abspath('pyCrossbill/database/pyCrossbill.sqlite')):
    # try to connect database 
    # 'pyCrossbill/database/pyCrossbill.sqlite' is the default database file on program
    try:
        conn = sqlite3.connect(database)
        cur = conn.cursor()
        # use sql command line to INSERT value or data from it.
        # '?' is use for load whatever data will be taken from data var. and collectivly call with data values
        query = """INSERT INTO NTTABLE (LOCUS, DEFINITION, VERSION, ORGANISM,
                REFRENCE, AUTHORS, COMMENTS, ORIGIN) VALUES (?,?,?,?,?,?,?,?);"""
        # load all the sql command from query var.
        data = (locus, definition, version, organism, refrence, authors, comments, origin)
        cur.execute(query, data) # execute above sql command 
        conn.commit()    # save all data to database table 
        cur.close()   # close cursor
        # return all succsess and database path
        return f'successfully {database}'
    # returns, if any condition of database failure
    except sqlite3.Error as err:        # sqlite3 error display
        return f'ERR: insertion Error {err}'
    

# all the program are same to above program 
# except this program is use for INSERT protein data to PTTABLE 
def InsertJob_prot(locus, definition, version, organism, refrence, authors, comments, origin, amino,
    database=os.path.abspath('pyCrossbill/database/pyCrossbill.sqlite')):
    try:
        conn = sqlite3.connect(database)
        cur = conn.cursor()
        query = """INSERT INTO PTTABLE (LOCUS, DEFINITION, VERSION, ORGANISM,
                REFRENCE, AUTHORS, COMMENTS, ORIGIN, AMINO_ACIDS) VALUES (?,?,?,?,?,?,?,?,?);"""
        # there is the another attribute amino is also use here that's defining number of amino acid molecule
        # presennt in sequence, all the procces are same to above program
        data = (locus, definition, version, organism, refrence, authors, comments, origin, amino)
        cur.execute(query, data)
        conn.commit()
        cur.close()
        return f'successfully {database}'
    except sqlite3.Error as err:
        return f'ERR: insertion Error {err}'


def search_DEFINITION_bySELECT(definition, table='protein', database=database):
    try:
        result = []
        con = sqlite3.connect(database)
        cur = con.cursor()
        if table == 'protein':
            query = ("""SELECT * FROM PTTABLE WHERE DEFINITION = '%s'""" %(definition))
        elif table == 'nucleotide':
            query = ("""SELECT * FROM NTTABLE WHERE DEFINITION = '%s'""" %(definition))
        
        try:
            cur.execute(query)
            for row in cur.fetchall():
                result.append(row)
            con.commit()
            cur.close()
            return result
        except UnboundLocalError:
            return "ERR: sorry, no any valid table is selected.. you can only accses 'protein' or 'nucleotide' table at a time"
        
    except sqlite3.Error as err:
        return f'searching Error occure: {err}'



########################################### END OF THE PROGRAM #####################################################

print(search_DEFINITION_bySELECT('this is virus', table='nucleotide'))
