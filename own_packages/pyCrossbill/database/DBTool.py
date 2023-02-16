#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 01:41:36 2021
database programs for storing fetching and analyse protein and nucleotide data
what's you are store
@author: ANIKET YADAV
"""

# By default pkg. use sqlite3 database 
import sqlite3
import os         # for defining path of the database file 

# specifying database ________
database = os.path.abspath('pyCrossbill/database/pyCrossbill.sqlite')

# insert data to nucleotide database table
# returns about the succses massege and path of the database file   
def InsertJob_nucleotide(locus, definition, version, organism, refrence, authors, comments, origin, 
    database = database):
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
    database=database):
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


# method of select use for searching data from table by definition 
# use default database for searching, and table you can select only 'protein' or 'nucleotide' at a same time
def search_DEFINITION_bySELECT(definition, table='protein', database=database):
    try:
        result = []
        con = sqlite3.connect(database)
        cur = con.cursor()
        # check if table protein is selected
        # so search or select from pre define PTTABLE by DEFINITION
        if table == 'protein':
            query = ("""SELECT * FROM PTTABLE WHERE DEFINITION = '%s'""" %(definition))
        # else table nucleotide is selected so search from NTTABLE by DEFINITION
        elif table == 'nucleotide':
            query = ("""SELECT * FROM NTTABLE WHERE DEFINITION = '%s'""" %(definition))
        
        # try to execute query and fetch all row from database
        try:
            cur.execute(query)
            for row in cur.fetchall():
                result.append(row)
            con.commit()
            cur.close()             # commit and close database after all the opration to be done ! 
            return result
        # return when you specify wrong table name or location
        except UnboundLocalError:
            return "ERR: sorry, no any valid table is selected.. you can only accses 'protein' or 'nucleotide' table at a time"
        
    except sqlite3.Error as err:
        return f'searching Error occure: {err}'


# program for printing whole data from a database
def ShowAll_database(table_name='nucleotide', database=database):
    try:
        result = []
        con = sqlite3.connect(database)
        cur = con.cursor()
        if table_name == 'protein':
            query = ("""SELECT * FROM PTTABLE""")
        elif table_name == 'nucleotide':
            query = ("""SELECT * FROM NTTABLE""")
        
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

# print(search_DEFINITION_bySELECT('this is virus', table='nucleotide'))
# print(InsertJob_nucleotide('28.38N', 'a bacterial genome', '26.4', 'bacteria', '2', 'aniket', 'nothing', 'AGTAGCGATCGGCTAGCTCGA'))
for i in ShowAll_database(table_name='nucleotide', database=database):
    print(i)
