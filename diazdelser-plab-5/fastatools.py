from termcolor import colored


def all_fasta_sequences(filename:str) -> list:
	"""Reads all sequences from fasta file and returns a list of tuples containing of header and sequence"""
	with open(filename, 'r') as f:
		# Read all lines
		all_records = "".join(f.readlines()).split('>')

		# Split into headers and sequences
		fasta_list = [ (record.split('\n',1)[0], record.split('\n',1)[1].replace('\n','')) \
					   for record in all_records if record !='']

		print(colored(f'Found a total of {len(fasta_list)} sequences in {filename}', 'green'))
	return fasta_list

def single_fasta_sequence(filename:str):
	"""Reads all sequences from fasta file and returns a list of tuples containing of header and sequence"""
	with open(filename, "r") as f:
		line = f.readline()
		while True:
			if line.startswith('>'):
				header = line.replace('\n','')
				# Read the rest of the lines as long as they're not headers
				seq = ''
				new_line = ''
				while not new_line.startswith('>'):
					seq = seq + str(new_line)
					try:
						new_line = next(f)
					except StopIteration:
						return
					line = new_line
				yield (header, seq.replace('\n',''))

def write_to_fasta(outfile:str, header:str, sequence:str):
	"""Writes the given sequence and its header to the output file in FASTA format"""
	print(f'>{header}', file=outfile)
	for i in range(0, len(sequence), 70):
		if sequence[i:i+70] != '':
			print(sequence[i:i+70], file=outfile)

def complementary(seq:str) -> str:
	"""Takes DNA sequence and returns complementary"""
	switch = { 'A' : 'T', 'T': 'A', 'C' : 'G', 'G': 'C'}
	# Use reverse to output it in 5' -> 3'
	return "".join(reversed([switch[nt] for nt in seq]))
