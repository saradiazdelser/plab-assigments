import click
import pandas as pd
from diazdelser_fastatools import get_sequence_positions


# ----------------------------
# FUNCTIONS
# ----------------------------


def delimit_size(pos_dict: dict, limit: int) -> dict:
	"""Remove positions that dont make the min size limit"""
	return {end: start for end, start in pos_dict.items() \
			if (int(end) - int(start) >= int(limit) or int(start) - int(end) >= int(limit))}


def correct_orfs(orf_fasta, genes_fasta, minimum):
	"""It should calculate and print comparsion info"""
	# Open and read each file

	with open(orf_fasta, 'r') as f:
		orf_pos = delimit_size(get_sequence_positions(f), minimum)
	with open(genes_fasta, 'r') as f:
		gen_pos = get_sequence_positions(f)

	# Compare:
	# Correctly predict stop
	match_stop = orf_pos.keys() & gen_pos.keys()
	# Correctly predict stop and start
	match = [orf_pos[x] for x in match_stop if orf_pos[x] == gen_pos[x]]
	# Non-coding ORFs
	mismatch = orf_pos.keys() - gen_pos.keys()

	if len(orf_pos) == 0:
		return {'n_orf': 0}

	return {'n_orf': len(orf_pos),
			'n_genes': len(gen_pos),
			'n_orf_correct': len(match),
			'ratio_orf_correct': len(match) / len(orf_pos),
			'n_stop_correct': len(match_stop),
			'ratio_stop_correct': (len(match_stop)/ len(orf_pos)),
			'n_incorrect': len(mismatch),
			'ratio_mismatch_correct': len(mismatch) / len(orf_pos),
			}


def comparison(orf_file, gene_file):
	"""Check results wiith different values for minimum gene length threshold"""
	list_of_dics = []
	for i in [50, 100, 150, 200, 250, 300, 350]:
		minimum = i * 3
		list_of_dics.append(correct_orfs(orf_file, gene_file, minimum))

	# Turn list of dics into dataframe
	df = pd.DataFrame(list_of_dics, index=[50, 100, 150, 200, 250, 300, 350])
	df.plot()


# ----------------------------
# MAIN
# ----------------------------

@click.group(invoke_without_command=True)
def cli():
	pass

@cli.command()
@click.argument('orf_file')
@click.argument('gene_file')
@click.argument('min')
def main(orf_file: str, gene_file: str, min: int):
	"""
	Calculates:
		 * total number of open reading frames
		 * total number of genes
		 * total number and ratio of open reading frames correctly predicting a gene
		 * total number and ratio of open reading frames correctly predicting at least the stop
		codon of a gene
		 * number of missed genes
	"""
	results = correct_orfs(orf_file, gene_file, min)
	df = pd.DataFrame(results.values(), index=results.keys(), columns=['Results'])
	print(df)
	return


@cli.command()
@click.argument('orf_file')
@click.argument('gene_file')
def compare(orf_file: str, gene_file: str):
	comparison(orf_file, gene_file)


if __name__ == '__main__':
	cli()
