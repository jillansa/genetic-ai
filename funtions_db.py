
"""
@author__ = "Juan Francisco Illan"
@license__ = "GPL"
@version__ = "1.0.1"
@email__ = "juanfrancisco.illan@gmail.com"
"""

import sqlite3

from classes import *

def initDb():
        
    con = sqlite3.connect('database/mydatabase.db')
    cursorObj = con.cursor()

    cursorObj.execute("CREATE TABLE execute ( id integer PRIMARY KEY AUTOINCREMENT, querry_seq text, creation_date date, k integer, match_score integer, mismatch_score integer, gap_score integer, seed_threshold integer, mode integer, numSecuencesDb integer, numSeedCalculate integer , numSecuencesCalculate integer, time_seed real, time_extends real)")
    con.close()

def deleteDb():

    con = sqlite3.connect('database/mydatabase.db')
    cursorObj = con.cursor()

    cursorObj.execute("DELETE from execute)")

    con.commit()
    con.close()

def store_data_execution(executeEntity):

    con = sqlite3.connect('database/mydatabase.db')
    cursorObj = con.cursor()
    cursorObj.execute('INSERT INTO execute (querry_seq, creation_date, k, match_score, mismatch_score, gap_score, seed_threshold, mode, numSecuencesDb, numSeedCalculate, numSecuencesCalculate, time_seed, time_extends) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', executeEntity)
    con.commit()
    con.close()

def get_statistic_execution(filterPlot):

    con = sqlite3.connect('database/mydatabase.db')
    cursorObj = con.cursor()
    sql = 'SELECT * FROM execute'

    if len(filterPlot) > 0:
        sql += ' WHERE ' + filterPlot

    cursorObj.execute(sql)
    rows = cursorObj.fetchall()
    
    listStatisticExecution = []

    for row in rows:
        se = StatisticExecution()
        se.id = row[0]
        se.querry_seq = row[1]
        se.creation_date = row[2]
        se.k = row[3]
        se.match_score = row[4]
        se.mismatch_score = row[5]
        se.gap_score = row[6]
        se.seed_threshold = row[7]
        se.mode = row[8] 
        se.numSecuencesDb = row[9]
        se.numSeedCalculate = row[10]
        se.numSecuencesCalculate = row[11]
        se.time_seed = row[12]
        se.time_extends = row[13]
        listStatisticExecution.append(se)

    con.close()
    return listStatisticExecution


