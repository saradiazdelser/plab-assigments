import numpy as np
from collections import Counter
from itertools import combinations
import pandas as pd
from termcolor import colored
import tabulate
import click

def read_file(file:str) -> list:
	"""Read file into list of lines"""
	with open (file, 'r') as f:
		return f.readlines()

def read_fasta(filename:str) -> str:
	"""Get sequence from fasta file into string"""
	with open(filename, 'r') as f:
		text = f.readlines()[1:]
	return "".join(text).replace("\n","")

def score_matrix_to_df(lines:list) -> (str, list):
	"""Read list of lines into dataframe"""
	lines = [line.strip() for line in lines]
	# Turn lines into list of scores
	lines =[ [ letter for letter in line.split() if letter !=' '] for line in lines ]
	# Separate header from rest of the strings
	header = lines[0]
	# Remove first item from lists (it's the aa)
	lines = [line[1:] for line in lines[1:]]
	# Generate pandas dataframe, remove emtpy rows
	df = pd.DataFrame(lines, columns=header).dropna(how='all')
	# Add aa symbols as row headers
	df = df.rename(dict(enumerate(header)))
	return df

def create_alignment_matrix(seqA:str, seqB:str,w, score) -> (np.ndarray,np.ndarray):
	"""Create a Neeedleman-Wunsch alignment matrix"""
	# Initialize: M(0; 0) = 0
	matrix = np.zeros((len(seqA),len(seqB)),dtype=int)
	traceback = np.zeros((len(seqA),len(seqB)),dtype=tuple)
	# first row M(0; j) = jw for j = 1 ... m,
	for i in range(len(seqA)-1):
		matrix[i,0] = i*w
	# first column M(i; 0) = iw for i = 1 ... n
	for j in range(len(seqB)-1):
		matrix[0,j] = j*w
		pass
	# Fill Rest of the matrix
	for i in range(len(seqA)-1):
		for j in range(len(seqB)-1):
			options = {
				matrix[i-1,j-1] + int(score.loc[seqA[i], seqB[j]]) : (i-1,j-1),
				matrix[i-1,j] + w : (i-1,j),
				matrix[i,j-1] : (i,j-1)
			}
			best_option = max(options.keys())
			matrix[i,j] = best_option
			traceback[i,j] = options.get(best_option)
	return matrix, traceback


def traceback(i:int,j:int,traceback_matrix:np.ndarray, seqA:str,seqB:str,alignement:list):
	"""Check traceback of Neeedleman-Wunsch alignment matrix"""
	if traceback_matrix[i,j] == (i-1,j-1):
		alignement.append((seqA[j],seqB[i]))
		traceback(i-1,j-1,traceback_matrix, seqA,seqB,alignement)
	elif traceback_matrix[i,j] == (i-1,j):
		alignement.append((seqA[j],'-'))
		traceback(i-1,j,traceback_matrix, seqA,seqB,alignement)
	elif traceback_matrix[i,j] ==(i,j-1):
		alignement.append(('-',seqB[i]))
		traceback(i,j-1,traceback_matrix, seqA,seqB,alignement)
	return alignement

def grouper(n, iterable):
	if iterable:
		if n < len(iterable):
			args = [iter(iterable)] * n
			return zip(*args)
		else:
			return [iterable]

def print_alignment(alignement:list,score_matrix:pd.DataFrame,n:int=80):
	"""Print sequence alignment"""
	# Print out a max of 80 characters per line
	zipper = grouper(n,iterable=alignement)
	if zipper:
		for alignement_bit in zipper:
			X,Y=zip(*alignement_bit)
			score_symbol = []
			for x,y in zip(X,Y):
				if x == y:
					score_symbol.append('|')
				elif (x=='-') or (y=='-'):
					score_symbol.append(" ")
				elif int(score_matrix.loc[x,y]) > 0:
					score_symbol.append(':')
				else:
					score_symbol.append(" ")

			print("".join(X))
			print("".join(score_symbol))
			print("".join(Y))
			print("\n")
	else:
		print("No alignment")


def total_score(alignement:list, score_matrix:pd.DataFrame) -> int or str:
	"""Calculate the total alignment score of a given alignment"""
	if alignement:
		x,y=zip(*alignement)
		scores = [ int(score_matrix.loc[x,y]) for x,y in zip(x,y) if (x != '-') and (y!= '-') ]
		return sum(scores)
	else:
		return 'N/A'


def Needleman_Wunsch(seqA:str, seqB:str, score_matrix:pd.DataFrame, w:int) -> int:
	"""Runs Neeedleman-Wunsch Algorithm"""
	# Construction of alignment matrix
	matrix, traceback_matrix = create_alignment_matrix(seqB,seqA, w,score_matrix)

	# Construct alignment via traceback matrix
	n = len(traceback_matrix)
	m = len(traceback_matrix[0])
	alignement = []
	alignement = traceback(n-2,m-2,traceback_matrix, seqA,seqB,alignement)
	print_alignment(alignement,score_matrix)
	return total_score(alignement, score_matrix)


# ----------------------------
# MAIN
# ----------------------------
@click.command()
@click.argument('fasta', nargs=-1)
@click.argument('score_file',nargs=1)
@click.argument('w',nargs=1)
def main(fasta:str, score_file:str,w:int=8):
	"""Use Needleman-Wunsch algorithm to create an alignement between the sequences in the inputed
	fasta files, using the inputed scoring matrix."""
	# Turn w into negative (click wont allow negative numbers as integers)
	w=-abs(int(w))
	# Open score matrix file
	text = read_file(score_file)
	# Read score into df
	score = score_matrix_to_df(text)

	# Open fasta files and read into list
	seqs = {}
	for file in fasta:
		seq = read_fasta(file)
		seqs.update({seq : file})

	# Calculate alignment for all combinations
	for seqA,seqB in combinations(seqs.keys(),2):
		print(colored(f'Alignment for: {seqs.get(seqA)} & {seqs.get(seqB)}\n', 'green'))
		# Try Needleman_Wunsch algorithm
		s= Needleman_Wunsch(seqA=seqA, seqB=seqB, score_matrix=score, w=w)
		print(colored(f'\nTotal score: {s}', 'green'))

if __name__ == '__main__':
	main()
