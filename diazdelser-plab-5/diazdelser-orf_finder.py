from termcolor import colored
import click
from diazdelser_fastatools import single_fasta_sequence, complementary, write_to_fasta, all_fasta_sequences
import re

# ----------------------------
# FUNTIONS
# ----------------------------

def orf_finder(seq:str) -> (list, int):
	# Start with ATG and end in TAA,TAG, or TGA.
	orf = r'(?:\w\w\w)*?((?:ATG)(?:\w\w\w)+?T(?:AG|AA|GA))'
	results = []
	sequences = []
	# Search every reading frame
	for i in range(0,3):
		for m in re.finditer(orf, seq[i:]):
			start, end = m.span(1)
			# Add matched group to results
			results.append((m.group(1), (start+i, end+i), end-start))
			sequences.append(m.group(1))
	return list(set(results))


# ----------------------------
# MAIN
# ----------------------------

@click.command()
@click.argument('inputfile')
@click.argument('outputfile')
def main(inputfile:str, outputfile:str):
	"""Find all the (longest) open reading frames in (from both strands of) the DNA given"""
	# Open input file and get sequence
	with open(inputfile, 'r') as infile:
		header, genome = all_fasta_sequences(infile)[0]
	# Open input file and add orf found
	with open(outputfile, 'w') as outfile:
		header = "|".join(header.split("|")[:4])
		# Length of seq
		length = len(genome)
		counter = 0
		# Find orf in forward strand
		for each in orf_finder(genome):
			# Add +1 to pos because python starts count at 0
			write_to_fasta(outfile, f'{header}|:{each[1][0]+1}-{each[1][1]}', each[0])
			counter+=1
		# Find orf in complementary strand
		for each in orf_finder(complementary(genome)):
			# Add +1 to pos because python starts count at 0
			write_to_fasta(outfile, f'{header}|:c{length-each[1][0]}-{length-each[1][1]+1}', each[0])
			counter+=1

	print(colored(f'Successfully saved {counter} ORFs to output file {outputfile}', 'green'))
	return

if __name__ == '__main__':
	main()