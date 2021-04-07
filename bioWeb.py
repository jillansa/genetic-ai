
"""
@author__ = "Juan Francisco Illan"
@license__ = "GPL"
@version__ = "1.0.1"
@email__ = "juanfrancisco.illan@gmail.com"
"""

from flask import Flask, request, render_template
from flask_debugtoolbar import DebugToolbarExtension
import sqlite3
import matplotlib.pyplot as plt
import random

from classes import *
from blast_process import *
from neuronal_process import *

app = Flask(__name__)
app.config['SECRET_KEY'] = "lailolailo"
app.debug = True
toolbar = DebugToolbarExtension(app)


# Controller to index
@app.route('/')
def home():
	return render_template('index.html')


# Controller to blast access
@app.route('/blast_nucleotide', methods=["GET"])
def blast_nucleotide():
	cb = ConfigBlast()
	return render_template('blast_nucleotide.html',configBlast=cb)

# Controller to blast request
@app.route('/blast_nucleotide', methods=["POST"])
def blast_nucleotide_post():

	# valores indicados en el formulario 
	cb = ConfigBlast()
	cb.querry_seq = request.form['querry_seq'] #The querry sequence we search for
	cb.k = int(request.form['k']) #Word length
	cb.match_score = int(request.form['match_score']) #Score added on match occurunce in alignment
	cb.mismatch_score = int(request.form['mismatch_score']) #Score added on mismatch occurunce in alignment
	cb.gap_score = int(request.form['gap_score']) #Score added on gap occurunce in alignment
	cb.seed_threshold = int(request.form['seed_threshold']) #The minimum score to considered a seed
	#cb.mode = int(request.form['mode']) # Mode of execution
	cb.num_secuences = int(request.form['num_secuences']) # Num of secuences to read from bd

	# read sequences in database
	db_sequences = read_data('database/db_seq.csv')
	# execute blast
	blastResult = blast_execute(cb,db_sequences[:cb.num_secuences+1])	
		
	#  consultamos ahora todas las stadisticas para mostrar un resumen de las ejeciones realizadas
	lse = get_statistic_execution('')    
	
	listSeedTimes = []
	listExtendsTimes = []
	
	for se in lse:
		listSeedTimes.append(se.time_seed)
		listExtendsTimes.append(se.time_extends)

	dataStats = {'seedTimes':listSeedTimes[len(listSeedTimes)-15:],
				'extendsTimes':listExtendsTimes[len(listExtendsTimes)-15:]}	
  
	hash = random.getrandbits(32)
	plt.figure()
	df = pd.DataFrame(dataStats)
	df.plot() 
	plt.xlabel('Nº of execution (last 15 execution)')
	plt.ylabel('Time of execution (sg.)')
	plt.savefig('static/images/temp/plot_' + str(hash) + '.png') 
	plt.close() 

	# response html with blastResult
	return render_template('blast_nucleotide.html', configBlast=cb, blastResult=blastResult, listStatisticExecution=lse, plotImage='plot_' + str(hash) + '.png')


# Controller to blast access
@app.route('/blast_protein', methods=["GET"])
def blast_protein():

	cbp = ConfigBlastProtein()	
	return render_template('blast_protein.html',configBlastProtein=cbp)


# Controller to blastProtein request
@app.route('/blast_protein', methods=["POST"])
def blast_protein_post():

	# valores indicados en el formulario 
	cbp = ConfigBlastProtein()
	cbp.querry_seq = request.form['querry_seq'] #The querry sequence we search for
	cbp.mode = int(request.form['mode']) # Mode of execution

	if (app.countVectorizer == None or app.classifierMNB == None) :
		seq_data =  pd.read_csv('database/db_seq.csv', sep=';', engine='python')	
		app.countVectorizer, app.classifierMNB, app.statsMNB = createClassifierMNB(seq_data)

	## mode process the secuences		
	blastResultProtein = clasiffierMNB(cbp, app.countVectorizer, app.classifierMNB)	

	clasiffierStats = app.classifierMNB
		
	# response html with blastResult
	return render_template('blast_protein.html', configBlastProtein=cbp, blastResultProtein=blastResultProtein, statsMNB=app.statsMNB)


# Controller to blast access
@app.route('/blast_pathogen', methods=["GET"])
def blast_pathogen():

	cbp = ConfigBlastPathogen()	
	return render_template('blast_pathogen.html',configBlastPathogen=cbp)


# Controller to blastProtein request
@app.route('/blast_pathogen', methods=["POST"])
def blast_pathogen_post():

	# valores indicados en el formulario 
	cbp = ConfigBlastPathogen()
	cbp.querry_seq = request.form['querry_seq'] #The querry sequence we search for
	cbp.mode = int(request.form['mode']) # Mode of execution

	if (app.tokenizerLSTM == None or app.modelLSTM == None) :
		seq_data =  pd.read_csv('database/pathogen_sars_cov_2.csv', sep=';', engine='python')		
		app.tokenizerLSTM, app.modelLSTM , app.statsLSTM, app.shortModelSummary = createClasiffierLSTM(seq_data)

	## mode process the secuences		
	blastResultPathogen = clasiffierLSTM(cbp, app.tokenizerLSTM, app.modelLSTM)	
		
	# response html with blastResult
	return render_template('blast_pathogen.html', configBlastPathogen=cbp, blastResultPathogen=blastResultPathogen, statsLSTM=app.statsLSTM, shortModelSummary = app.shortModelSummary)

# Controller to blast request
@app.route('/openAxisPlot', methods=["POST"])
def openAxisPlot():
	axisPlot = request.form['AxisPlot'] # The mode plot to generate
	filterPlot = request.form['filterPlot'] # The filter for plot generate
	
	#  consultamos ahora todas las stadisticas para mostrar un resumen de las ejeciones realizadas
	lse = get_statistic_execution(filterPlot)   
	listSeedTimes = []
	listExtendsTimes = []
	listLenQuery = []
	listKValue = []
	listNumSequences = []
	listMode = []

	for se in lse:
		listSeedTimes.append(se.time_seed)
		listExtendsTimes.append(se.time_extends)
		listNumSequences.append(se.numSecuencesDb)
		listLenQuery.append(len(se.querry_seq))
		listKValue.append(se.k)
		listMode.append(se.mode)

	if int(axisPlot)==1: # 1 - creation_date
		dataStats = {'seedTimes':listSeedTimes,
					'extendsTimes':listExtendsTimes}	
						
		df = pd.DataFrame(dataStats)
		df.plot() # x ='time_seed', y='time_extends', kind = 'scatter'
		plt.xlabel('Nº of execution')
		plt.ylabel('Time of execution (sg.)')

	elif int(axisPlot)==2: # 2 - length(querry_seq) 
		f = plt.figure()    
		f, axes = plt.subplots(nrows = 1, ncols = 2)
		
		axes[0].scatter(listLenQuery, listSeedTimes, marker = "x", color='blue')
		axes[0].set_xlabel('Length(query) of execution')
		axes[0].set_ylabel('Time of execution SeedTimes (sg.)')

		axes[1].scatter(listLenQuery, listExtendsTimes, marker = 'x', color='red')
		axes[1].set_xlabel('Length(query) of execution')
		axes[1].set_ylabel('Time of execution ExtendsTimes (sg.)')
		axes[1].yaxis.set_label_position("right")

	elif int(axisPlot)==3: # 3 - k
		f = plt.figure()    
		f, axes = plt.subplots(nrows = 1, ncols = 2)
		
		axes[0].scatter(listKValue, listSeedTimes, marker = "x", color='blue')
		axes[0].set_xlabel('K-parameter of execution')
		axes[0].set_ylabel('Time of execution SeedTimes (sg.)')

		axes[1].scatter(listKValue, listExtendsTimes, marker = 'x', color='red')
		axes[1].set_xlabel('K-parameter of execution')
		axes[1].set_ylabel('Time of execution ExtendsTimes (sg.)')
		axes[1].yaxis.set_label_position("right")
	elif int(axisPlot)==4: # 4 - mode
		f = plt.figure()    
		f, axes = plt.subplots(nrows = 1, ncols = 2)
		
		axes[0].scatter(listMode, listSeedTimes, marker = "x", color='blue')
		axes[0].set_xlabel('Mode of execution')
		axes[0].set_ylabel('Time of execution SeedTimes (sg.)')

		axes[1].scatter(listMode, listExtendsTimes, marker = 'x', color='red')
		axes[1].set_xlabel('Mode of execution')
		axes[1].set_ylabel('Time of execution ExtendsTimes (sg.)')
		axes[1].yaxis.set_label_position("right")

	elif int(axisPlot)==5: # 5 - numSecuencesDb
		f = plt.figure()    
		f, axes = plt.subplots(nrows = 1, ncols = 2)
		
		axes[0].scatter(listNumSequences, listSeedTimes, marker = "x", color='blue')
		axes[0].set_xlabel('NumSequences in db')
		axes[0].set_ylabel('Time of execution SeedTimes (sg.)')

		axes[1].scatter(listNumSequences, listExtendsTimes, marker = 'x', color='red')
		axes[1].set_xlabel('NumSequences in db')
		axes[1].set_ylabel('Time of execution ExtendsTimes (sg.)')
		axes[1].yaxis.set_label_position("right")
	
	hash = random.getrandbits(32)
	plt.savefig('static/images/temp/plot_' + str(hash) + '.png') 
	plt.close()  

	# response html with plotimage
	return render_template('plot.html',plotImage='plot_' + str(hash) + '.png', filterPlot=filterPlot)

# Controller to blast access
@app.route('/about', methods=["GET"])
def about():

	return render_template('about.html')

# startup HTTP web service
if __name__ == '__main__':
	print("Running service ...")
	app.classifierMNB = None
	app.countVectorizer = None
	app.tokenizerLSTM = None
	app.modelLSTM = None
	app.run(host='0.0.0.0')

