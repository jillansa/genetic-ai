"""
@author__ = "Juan Francisco Illan"
@license__ = "GPL"
@version__ = "1.0.1"
@email__ = "juanfrancisco.illan@gmail.com"
"""

class ConfigBlast:
	querry_seq = ''
	k = 0
	match_score = 0
	mismatch_score = 0
	gap_score = 0
	seed_threshold = 0
	mode = 0
	num_secuences = 0
	
	def __init__(self):
		self.querry_seq = 'AAATCTGTTCGCTTCATT' #The querry sequence we search for
		self.k = 5 #Word length
		self.match_score = 2 #Score added on match occurunce in alignment
		self.mismatch_score = -1 #Score added on mismatch occurunce in alignment
		self.gap_score = -1 #Score added on gap occurunce in alignment
		self.seed_threshold = 6 #The minimum score to considered a seed
		self.num_secuences = 3 # Num of secuences to read from bd

class ConfigBlastProtein:
	querry_seq = ''

	def __init__(self):
		self.querry_seq = 'ACGGTGCTACTCAAGGCCCGGGAAGGTGGCGGTGGAAATCGCAAAGGCAAAAGCAAGAAATGGCGGCAGATGC' #The querry sequence we search for


class ConfigBlastPathogen:
	patogen = ''
	querry_seq = ''

	def __init__(self):
		self.patogen = 'SARS CoV 2'
		self.querry_seq = 'ACGGCAGTGAGGACGATCAGACAACTACTATTCAAACAATGTTGAGGTTCAACCTCAATTAGAGATGGAACTTACACAGTTGTTCAGACTATTGAAGTGAATAGTTTTA' #The querry sequence we search for


class BlastResult:
	secuences = []	# List of AlignmentSecuence
	time_seed = None
	time_extends = None

class AlignmentSecuence:
	idSecuence = None 
	nameSecuence = None 
	strSecuence = None 
	scoreMaxSecuence = None 
	alignments = []	# List of Alignments calculate for Secuence

	def __lt__(self, other):
		if self.scoreMaxSecuence < other.scoreMaxSecuence:
			return True
		else:
			return False

class StatisticExecution:
	id = None
	querry_seq = None
	creation_date = None
	k = None
	match_score = None
	mismatch_score = None
	gap_score = None
	seed_threshold = None
	mode = None
	numSecuencesDb = None
	numSeedCalculate = None
	numSecuencesCalculate = None
	time_seed = None
	time_extends = None
