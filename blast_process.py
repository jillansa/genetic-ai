
"""
@author__ = "Juan Francisco Illan"
@license__ = "GPL"
@version__ = "1.0.1"
@email__ = "juanfrancisco.illan@gmail.com"
"""

import pandas as pd
import numpy as np
import math
import re
import time
import sqlite3
from datetime import date, datetime

from classes import *
from funtions_db import *


""" Read the database in filepath and return a list of sequences found in database """
def read_data(filepath):
    reader = open(filepath,'r')
    return re.split('\n',reader.read()[1:])

    
""" The algorithm for a request of alignment query over database secuences"""    
def blast_execute(cb, db_sequences):

    querry_seq = cb.querry_seq
    k = cb.k
    match_score = cb.match_score
    mismatch_score = cb.mismatch_score
    gap_score = cb.gap_score
    seed_threshold = cb.seed_threshold

    blastResult = BlastResult()
    blastResult.secuences.clear()

    time_seed = 0.0
    time_extends  = 0.0
    num_seed_alignment = 0
    num_seq_alignment = 0

    # For each secuence in db 
    for i in range(len(db_sequences)):

        # skip first row db_secuence
        if i == 0:
            continue

        id_seq = []
        names = []
        querry_seed = []
        query_index_seed = []
        db_seed = []
        db_index_seed = []
        score_seed = []
        querry_alignments  = []
        db_alignments  = []
        row_scores  = []
        row_seed_max_score = []
    
        db_seq = db_sequences[i].split(';')[6]
        
        # find the seed over db_seq 
        start_seed_time = time.time()
        #if (cb.mode == 2):
        seed_table = find_seeds(querry_seq,db_seq,k,match_score,mismatch_score,seed_threshold)
        #else:
        #seed_table = find_seeds_trie(querry_seq,db_seq,k,match_score,mismatch_score,seed_threshold)

        final_seed_time = time.time()
        time_seed = time_seed + final_seed_time - start_seed_time

        num_seed_alignment = num_seed_alignment + seed_table.shape[0]

        # extends seeds with gaps 
        start_extends_time = time.time()
        table_extend = extend_seed(seed_table,querry_seq,db_seq,k,match_score,mismatch_score,gap_score)
        final_extends_time = time.time()
        time_extends = time_extends + final_extends_time - start_extends_time
    
        score_max_extend_seed = 0
        j=0
        
        # si hemos obtenido alineamientos como resultado para la secuencia db_seq
        if len(table_extend.index) > 0:
                    
            for j in range(len(table_extend.index)):
            
                # guardamos los datos en la listas
                id_seq.append(db_sequences[i].split(';')[4])
                names.append(db_sequences[i].split(';')[5])
                querry_seed.append(table_extend.at[j, 'querry_seed'])
                query_index_seed.append(table_extend.at[j, 'query_index_seed'])
                db_seed.append(table_extend.at[j, 'db_seed'])
                db_index_seed.append(table_extend.at[j, 'db_index_seed'])
                score_seed.append(table_extend.at[j,'score_seed'])
                querry_alignments.append(table_extend.at[j, 'querry_alignment_extends'])
                db_alignments.append(table_extend.at[j, 'db_alignment_extends'])
                row_scores.append(table_extend.at[j, 'row_score'])

                if (table_extend.at[j, 'row_score'] > score_max_extend_seed):
                    score_max_extend_seed = table_extend.at[j, 'row_score']
                
                num_seq_alignment = num_seq_alignment + 1
            
            if (j>0):
                for q in range(j+1):
                    row_seed_max_score.append(score_max_extend_seed)

            # Creamos el data con los array de todos los alineamientos calculados para esa semilla
            dataAlignment = {
                #'id_seq':id_seq,
                #'names':names,
                'querry_seed':querry_seed,
                'query_index_seed':query_index_seed,
                'db_seed':db_seed,
                'db_index_seed':db_index_seed,
                'score_seed':score_seed,
                'querry_alignment_extends':querry_alignments,
                'db_alignment_extends':db_alignments,        
                'row_scores':row_scores,
                #'row_seed_max_score':row_seed_max_score
                }

            alignmentSecuence = AlignmentSecuence()
            alignmentSecuence.alignments.clear()

            alignmentSecuence.idSecuence = db_sequences[i].split(';')[4]
            alignmentSecuence.nameSecuence = db_sequences[i].split(';')[5]
            alignmentSecuence.strSecuence = db_seq
            alignmentSecuence.scoreMaxSecuence = score_max_extend_seed
            
            df = pd.DataFrame(dataAlignment)
            dfSort = df.sort_values(by=['row_scores'],ascending = [False]).head(100)
            alignmentSecuence.alignments = dfSort.values.tolist()           
            
            blastResult.secuences.append(alignmentSecuence)        
    
    # una vez finalizado el algoritmo, ordenamos para mostrar las mejores secuencias primero    
    blastResult.secuences.sort(reverse=True)
    
    # tomamos muestras del tiempo en la ejecucion de cada etapa
    blastResult.time_seed = time_seed
    blastResult.time_extends = time_extends

    # guardamos las estadistica en bd 
    executeEntity = (querry_seq, date.today(), k, match_score, mismatch_score, gap_score, seed_threshold, 1, len(db_sequences), num_seed_alignment, num_seq_alignment, blastResult.time_seed, blastResult.time_extends)
    store_data_execution(executeEntity)

    return blastResult


""" Finds the seeds with score > seed_threshold """
def find_seeds(querry_seq,db_seq,k,match_score,mismatch_score,seed_threshold):   
    
    # TODO mejora propuesta: buscar las semillas tras un preprocesamiento inicial que 
    # permita la busqueda agil de todas las cadenas que contiene conincidencia con una semilla
    # actuando como una estructura de acceso rapido hash, arbol trie, ...

    # split querry_seq in kmers 
    querry_kmers = extract_kmers(querry_seq,k)

    # split db_seq in kmers    
    db_kmers = extract_kmers(db_seq,k)

    kmers =[]
    querry_kmer =[]
    q_indicies=[]
    db_indicies=[]
    scores=[]
    
    # for each kmers in query_seq
    for i in range(len(querry_kmers)):
        # for each kmers in db_seq
        for j in range(len(db_kmers)):
            # calculate score
            score = ungapped_alignment(querry_kmers[i],db_kmers[j],match_score,mismatch_score)
            # if > seed_threshold, is a valid seed 
            if(score>=seed_threshold):
                kmers.append(db_kmers[j])
                q_indicies.append(i) 
                db_indicies.append(j)
                scores.append(score)
                querry_kmer.append(querry_kmers[i])
    
    data = {'db_kmer':kmers,
            'querry_index':q_indicies,
            'db_index':db_indicies,
            'score': scores,
            'querry_kmer': querry_kmer}
    seed_table = pd.DataFrame(data,columns = ["db_kmer","querry_index","db_index","score", "querry_kmer"])        
    return seed_table

    
""" Extracts every possible kmer from a sequence."""
def extract_kmers(sequence,k):
    kmers=[]
    kmer = ""
    for i in range(len(sequence)-k+1):
        kmer = sequence[i]
        for j in range(1,k):
            kmer+=sequence[i+j]      
        kmers.append(kmer)
    return kmers

""" Calculate ungapped alignment between two kmers. """
def ungapped_alignment(kmer1,kmer2,match_score,mismatch_score):
    scores = np.zeros(len(kmer1))
    for i in range(len(kmer1)):
        if(kmer1[i]==kmer2[i]): 
            scores[i] = scores[i-1]+match_score
        else:
            scores[i] = scores[i-1]+mismatch_score

    return scores[-1]


""" Extends the seeded with gaps """
def extend_seed(seed_table,querry_seq,db_seq,k,
                match_score,mismatch_score,gap_score):
    querry_seed = []
    querry_alignments = []
    db_seed = []
    db_alignments = []
    row_scores = []
    score_seed=[]
    query_index_seed=[]
    db_index_seed=[]

    # for each seed detected
    for i in range(seed_table.shape[0]): # seed_table.shape[0] : number of seed 
        # TODO en lugar de utilizar iloc, utilizar un acceso a traves de clave
        s,q,d = SmithWatermanMatrix(querry_seq,db_seq,seed_table.iloc[i,1],seed_table.iloc[i,2],k,
                    match_score,mismatch_score,gap_score,
                    seed_table.iloc[i,3])
        querry_alignments.append(q)        
        db_alignments.append(d)
        row_scores.append(s)

        # data of original seed to debug algoritm 
        db_seed.append(seed_table.iat[i,0]) # iat[row,column=0] db_kmer
        query_index_seed.append(seed_table.iat[i,1]) #query_index
        db_index_seed.append(seed_table.iat[i,2]) #db_index     
        score_seed.append(seed_table.iat[i,3]) #score_seed
        querry_seed.append(seed_table.iat[i,4]) #query_seed
        
    data = {'querry_seed':querry_seed,
            'query_index_seed':query_index_seed,
            'db_seed':db_seed,
            'db_index_seed':db_index_seed,
            'score_seed':score_seed,
            'querry_alignment_extends':querry_alignments,
            'db_alignment_extends':db_alignments,
            'row_score':row_scores}
    df = pd.DataFrame(data)    
    return df 


""" Process of extends with Smith-Waterman algorithm."""
def SmithWatermanMatrix(querry_seq,db_seq,querry_index,db_index,k,
        match_score,mismatch_score,gap_score,
        initial_score):
   
   querry_alignment = ""
   db_alignment = ""
   row_score = 0
   max_score = initial_score
   
   # create matrix SWM
   matrix = np.zeros((len(querry_seq),len(db_seq)))
   # initial value for seed evaluate to extend
   matrix[querry_index+k-1][db_index+k-1] = initial_score
   if (len(querry_seq) > querry_index+k):
      matrix[querry_index+k][db_index+k-1] = initial_score + gap_score
   if (len(db_seq) > db_index+k):
      matrix[querry_index+k-1][db_index+k] = initial_score + gap_score
   row = 0
   col = 0

   cont = True
   
   # Start building table from the kmer position
   # Inicializate index i,j
   i = querry_index+k
   j = db_index+k

   while (i < len(querry_seq)):
        #if(not cont):
        #   cont = True
        #   break

        # evitar que nos salgamos del recorrido en i sobre querry seq
        if (i >= len(querry_seq)):   
           break
        
        j_aux = j
        while (j < len(db_seq)):
           #if(not cont):
           #    cont = True
           #    break
           R = 0
           C = 0

           if(querry_seq[i] == db_seq[j]):
               D = matrix[i-1][j-1] + match_score
           else:
               # Evaluate mismatch
               D = matrix[i-1][j-1] + mismatch_score

           # Evaluate Gap in db sequence
           R = matrix[i][j-1] + gap_score
           # Evaluate Gap in query sequence
           C = matrix[i-1][j] + gap_score
               
           matrix[i][j] = np.max([D,R,C])
           row_score = matrix[i][j] 
           
           # Mejor puntuacion hasta el momento para esta seed
           if (row_score > max_score):
               max_score = row_score
               row = i
               col = j
           
           # Acotar recoridos de la matriz, si hay match y es solucion optima, avanzar en la diagonal
           if(querry_seq[i] == db_seq[j]):
               if (D == np.max([D,R,C])):
                   if (i+1 < len(querry_seq)):
                       matrix[i+1][j] = matrix[i][j] + gap_score
                   if (j+1 < len(db_seq)):
                       matrix[i][j+1] = matrix[i][j] + gap_score
                   j+=1
                   break
            
           # Si no hay match, puede llegar el punto que no tenga sentido seguir evaluando posibles GAPs
           if(querry_seq[i] != db_seq[j] and j > j_aux):  
               if (matrix[i][j] + match_score < max_score):
                j = j_aux
                break
           
           j+=1

           # j, si cae muy por debajo del score inicial, no tiene sentido seguir valorando esa fila sobre j++   
           #Break if we drop below threshold
           #if(i>= start_checking):
           #if(row_score < valorLimite):
               #continue = False
               #row = i
               #col = j
               #break
           #else:
           #= row_score 

        # TODO acortar la i en cada fila por delante u por detras
        #if(row_score < valorLimite):
           #continue = False
           #row = i
           #col = j
           #break
        i+=1
   # finalizada la extension de la semilla
   
   # reconstruimos desde la mejor puntuacion matrix[row][col] hacia atras  
   while(row >= querry_index+k and col>=db_index+k):
        # retroceder match
        if(matrix[row][col] == matrix[row-1][col-1]+match_score):
            querry_alignment += querry_seq[row]
            db_alignment += db_seq[col]
            row-=1
            col-=1
        # retroceder mistmach
        elif(matrix[row][col] == matrix[row-1][col-1]+mismatch_score):
            querry_alignment += querry_seq[row].lower()
            db_alignment += db_seq[col].lower()
            row-=1
            col-=1   
        # retroceder GAP de filas
        elif(matrix[row][col] == matrix[row-1][col]+gap_score):
            querry_alignment += querry_seq[row].lower()
            db_alignment += "_"
            row -=1
        # retroceder GAP de columnas
        elif(matrix[row][col] == matrix[row][col-1]+gap_score):
            db_alignment += db_seq[col].lower()
            querry_alignment += "_"
            col -=1  
        else:
            querry_alignment += 'X'
            db_alignment += 'X'
            row-=1
            col-=1    

   # invertimos las cadena
   querry_alignment = querry_alignment[::-1]
   db_alignment = db_alignment[::-1]                
   
   querry_seq_aux=''
   db_seq_aux=''
   # añadimos a las cadenas la parte de la semilla 
   for index in range(k):
      if(querry_seq[querry_index+index]==db_seq[db_index+index]): 
         querry_seq_aux += querry_seq[querry_index+index]
         db_seq_aux += db_seq[db_index+index]
      else:
         querry_seq_aux += querry_seq[querry_index+index].lower()
         db_seq_aux += db_seq[db_index+index].lower()

   querry_alignment = querry_seq_aux + querry_alignment
   db_alignment = db_seq_aux + db_alignment
    
   return max_score,querry_alignment,db_alignment

