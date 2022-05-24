import numpy as np
from collections import Counter
from itertools import combinations
import pandas as pd
from termcolor import colored
import click

# --------------------------------
# STATS FUNCTIONS
# --------------------------------

# 20 aa
aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def fa(seq:str) -> dict:
	"""Return frequency of all amino acids"""
	# Start count at 1
	counter = Counter(seq)
	add_1 = {key :1 for key in list(counter.keys()) }
	return counter + Counter(add_1)

def pair_freq(seqA,seqB):
	"""Add pair frequencies from two sequences into dictionary"""
	all_pairs = []
	# seq1 should be the smaller one
	if len(seqA) == len(seqB):
		seq1, seq2 = seqA, seqB
	else:
		seq1 = min([seqA,seqB],key = len)
		seq2 = max([seqA,seqB], key = len)

	for i in range(len(seq1)):
		# Append pair (ordered) to list:
		key = sorted([seq1[i],seq2[i]])
		all_pairs.append("".join(key))
	return all_pairs

def all_pair_freq(seq_list:list) -> dict:
	"""Get pair frequencies from list of sequences"""
	all_pairs = []
	for seq1,seq2 in combinations(seq_list,2):
		all_pairs.extend(pair_freq(seq1,seq2))

	counter = Counter(all_pairs)
	add_1 = {key :1 for key in list(counter.keys()) }
	return counter + Counter(add_1)

def rel_freq(frequency:dict) -> dict:
	"""Determine relative frequency for each aminoacid"""
	return { key : (val/sum(frequency.values())) for key, val in frequency.items() }

def expected_prob(pab:dict) -> dict:
	"""Determine the expected probability of each aminoacid"""
	expected = {}
	# Get all possible pairs
	possible_pairs = ["".join(sorted([aa1,aa2])) for aa1,aa2 in combinations(pab.keys(),2) ]
	# add pairs of the same aa
	possible_pairs.extend([aa+aa for aa in pab.keys() ])
	for pair in possible_pairs:
		# eaa = pa * pa
		if pair[0]==pair[1]:
			prob = pab.get(pair[0],0)*pab.get(pair[0],0)
		# eab = pa * pb + pb * pa = 2 * pa * pb
		else:
			prob = 2*pab.get(pair[0],0)*pab.get(pair[1],0)
		expected.update({pair: prob})
	return expected

def pair_score(pab, eab) -> dict:
	"""Calculate the score for each aminoacid pair"""
	scores = {}
	for pair in pab.keys():
		# sab = 2 * log2(pab/eab)
		score = 2*np.log2(pab.get(pair,1)/eab.get(pair,1))
		scores.update({pair : round(score,0)})
	return scores

# --------------------------------
# SCRIPT FUNCTIONS
# --------------------------------

# (a) read in the alignment data in an appropriate data structure
def read_file(file:str) -> list:
	"""Read file into list of lines"""
	with open (file, 'r') as f:
		return f.readlines()

def pre_processing(seqs:list) -> list:
	"""Pre-process alignment data:
	- Add padding at the end of shorter strings to make sure they're the same size
	- Remove newlines
	"""
	max_length = max([len(each) for each in seqs ])
	seqs = [ seq.ljust(max_length).replace('\n','') for seq in seqs ]
	return seqs


# (b) determine the log-odds scores for each possible alignment of amino acids
def calc_score(seqs:list) -> dict:
	"""Determine score from list of pre-processed aligned sequences"""

	freq = fa("".join(seqs))
	f_ab = all_pair_freq(seqs)

	p_a = rel_freq(freq)
	p_ab = rel_freq(f_ab)

	e_ab = expected_prob(p_a)
	p_score = pair_score(p_ab, e_ab)

	return p_score

# (c) produce a (nicely formatted) output of the resulting scoring matrix.
def calc_score_matrix(score:dict,aa) -> np.ndarray:
	"""Generate a score matrix from a score dictionary"""
	# Create matrix of 0 of size (20,20)
	matrix = np.zeros((20,20),dtype=int)

	# Assign index to each aa
	aa_index = { key: i for i,key in enumerate(aa)}
	# Iterate through score dictionary and input values into matrix
	for key, val in score.items():
		matrix[aa_index[key[0]],aa_index[key[1]]] = val
		matrix[aa_index[key[1]],aa_index[key[0]]] = val

	return matrix

def score_table(score_matrix:np.ndarray,aa, file) -> pd.DataFrame:
	"""Turn a numpy score matrix into a dataframe"""
	df = pd.DataFrame(score_matrix, index=aa, columns=aa)

	with open(file, 'w') as f:
		print(df.to_string(), file=f)
		print(colored(f'Score matrix saved to "{file}"'), 'green')

	return df


# --------------------------------
# MAIN
# --------------------------------

@click.command()
@click.argument('alignment')
@click.argument('output')
def main(alignment:str, output:str):
	"""Generates a score matrix in a text file."""
	seqs = read_file(alignment)
	seqs = pre_processing(seqs)
	score = calc_score(seqs)
	score_matrix = calc_score_matrix(score, aa)
	score_table(score_matrix, aa, output)

if __name__ == '__main__':
	main()
