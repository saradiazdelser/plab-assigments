import click
import pandas as pd
from diazdelser_fastatools import get_sequence_positions


# ----------------------------
# FUNCTIONS
# ----------------------------


def delimit_size(pos_dict: dict, limit: int) -> tuple[dict, dict]:
	"""Remove positions that dont make the min size limit"""
	big, small = {}, {}
	for end, start in pos_dict.items():
		if (int(end) - int(start) >= int(limit) or int(start) - int(end) >= int(limit)):
			big[end] = start
		else:
			small[end] = start

	return big, small

def compare_orfs(orf_pos, gen_pos) -> tuple:
	"""Find full matches (STOP & START), half-matches (STOP) and
	mismatches for positions in ORF and gene dictionaries"""
	# Correctly predict stop
	match_stop = orf_pos.keys() & gen_pos.keys()
	# Correctly predict stop and start
	match = [orf_pos[x] for x in match_stop if orf_pos[x] == gen_pos[x]]
	# Non-coding ORFs
	mismatch = orf_pos.keys() - gen_pos.keys()
	return match, mismatch, match_stop


def calculate_acc(TP:int, FN:int, TN:int, FP:int) -> float:
	"""0,5 * ( TP/(TP+FN) + TN/(TN+FP) )"""
	return (0.5 * ( TP/(TP+FN) + TN/(TN+FP) ))

def correct_orfs(orf_fasta, genes_fasta, minimum):
	"""It should calculate and print comparsion info"""
	# Open and read each file

	with open(orf_fasta, 'r') as f:
		orf_pos, orf_pos_small = delimit_size(get_sequence_positions(f), minimum)
	with open(genes_fasta, 'r') as f:
		gen_pos = get_sequence_positions(f)

	# Compare:
	match, mismatch, match_stop = compare_orfs(orf_pos, gen_pos)
	match_small, mismatch_small, match_stop_small = compare_orfs(orf_pos_small, gen_pos)

	if len(orf_pos) == 0 or len(orf_pos_small) ==0:
		return {'n_orf': 0}

	acc = {'TP':len(match),
		   'FP':len(mismatch)+len(match_stop),
		   'TN':len(mismatch_small)+len(match_stop_small),
		   'FN':len(match_small)
		}
	acc_balance = calculate_acc(acc['TP'], acc['FN'], acc['TN'], acc['FP'])


	results = {'n_orf': len(orf_pos),
			'n_genes': len(gen_pos),
			'n_orf_correct': len(match),
			'ratio_orf_correct': len(match) / len(orf_pos),
			'n_stop_correct': len(match_stop) - len(match),
			'ratio_stop_correct': (len(match_stop) - len(match)) / len(orf_pos),
			'n_incorrect': len(mismatch),
			'ratio_mismatch_correct': len(mismatch) / len(orf_pos),
			'acc': acc_balance}
	return results


def comparison(orf_file, gene_file):
	"""Check results wiith different values for minimum gene length threshold"""
	list_of_dics = []
	for i in [50, 100, 150, 200, 250, 300, 350]:
		minimum = i * 3
		results = correct_orfs(orf_file, gene_file, minimum)
		list_of_dics.append(results)


	# Turn list of dics into dataframe
	df = pd.DataFrame(list_of_dics, index=[50, 100, 150, 200, 250, 300, 350])
	return df


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
@click.option('--acc', is_flag=True)
def compare(orf_file: str, gene_file: str, acc:bool):
	df = comparison(orf_file, gene_file)
	if acc:
		df['acc'].plot()
	else:
		df.drop(['acc'], axis=1).plot()


if __name__ == '__main__':
	cli()
