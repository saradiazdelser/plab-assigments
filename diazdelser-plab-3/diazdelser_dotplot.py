import numpy as np
import matplotlib.pyplot as plt
import click
from pathlib import Path
from termcolor import colored

# (4 pts) Generate dotplot matrix
def dotplot(seqA:str,seqB:str,w:int,s:int) -> np.ndarray:
    """Returns a dotplot matrix of matches between two input sequences"""
    dp = np.zeros((len(seqA), len(seqB)),dtype=int)
    k = int((w-1)/2)
    # Add spaces too the begining of both sentences to slice strings
    seqA = seqA.rjust(len(seqA)+k, ' ')
    seqB = seqB.rjust(len(seqB)+k, ' ')
    # Parse through both sequences
    for i in range(len(seqA)):
        for j in range(len(seqB)):
            if seqA[i] == seqB[j] and seqA[i] != ' ':
                # Check if window has s matches
                score = 0
                for x, y in zip(seqA[i-k:i+k], seqB[j-k:j+k]):
                    if x == y:
                        score+=1
                if score >= s:
                    # Assign 1 to spot
                    dp[i-k,j-k] = 1

    return dp

# (4 pts) The Dotplot as ASCII art
def dotplot2Ascii(dp:np.ndarray, seqA:str, seqB:str,heading:str,filename:str):
    """Creates an output file for the doplot"""
    with open (filename, 'w') as file:
        print(heading,"\n", file=file)
        print(f' |{seqA}', file=file)
        print(f'-+{"-"*len(seqA)}', file=file)
        for j in range(len(dp[0])):
            row = "".join([str(each) for each in dp[:,j]]).replace('1','*').replace('0',' ')
            print(f'{seqB[j]}|{row}', file=file)
    return


# (4 pts) Graphical output using matplotlib
def dotplot2Graphics(dp:np.ndarray,labelA:str,labelB:str,heading:str,filename:str,spFlag:bool=False,bigFlag:bool=False):
    """Creates a dotplot, saves it to a file and displays it on screen."""
    # -----------------
    # Type of graph
    # -----------------

    if spFlag:
        # Get X and Y coordinates from dotplot matrix
        x,y = [],[]
        for i in range(len(dp)):
            for j in range(len(dp[0])):
                if dp[i][j] == 1:
                    x.append(i)
                    y.append(j)
        # Create figure with scatter
        fig = plt.scatter(x,y, cmap="Blues", marker="+", s=1)

    else:
        # Create figure with imshow
        fig = plt.imshow(dp, cmap='Blues')


    # -----------------
    # Size of graph
    # -----------------
    # if the image is big, remove the letters from the sides of the dp
    if bigFlag:
        # Add labels
        plt.xlabel(labelA)
        plt.ylabel(labelB)

        # Hide ticks, show only letters
        fig.axes.tick_params(top=False, bottom=False, left=False, right=False)
        fig.axes.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    # if the image is small, add each letter to the corresponding tick
    else:
        # Add letters to ticks
        plt.xticks(np.arange(0,len(labelA),1))
        plt.yticks(np.arange(0,len(labelB),1))
        fig.axes.set_xticklabels(list(labelA))
        fig.axes.set_yticklabels(list(labelB))

        # Hide ticks, show only letters
        fig.axes.tick_params(top=False, bottom=False, left=False, right=False)
        fig.axes.tick_params(labeltop=True, labelbottom=False)

    plt.axis("tight")  # gets rid of white border
    plt.axis("image")  # square up the image instead of filling the "figure" space


    # Add title
    plt.title(heading,fontweight="bold",pad=10)
    plt.savefig(filename)

    plt.show()
    return


# (4 pts) Write a script
def read_fasta(filename:str) -> str:
    """Get sequence from fasta file into string"""
    with open(filename, 'r') as f:
        text = f.readlines()[1:]
    return "".join(text).replace("\n","")



@click.command()
@click.argument('w')
@click.argument('s')
@click.argument('seqA')
@click.argument('seqB')
@click.argument('title')
@click.argument('output')
def main(w:int, s:int, seqa:str, seqb:str, title:str, output:str):
    """Generates a sequence match dotplot."""
    # Read fasta
    print(colored(f'Reading files: \n\t{seqa}\n\t{seqb}', 'yellow'))
    seqA_string = read_fasta(filename=seqa)
    seqB_string = read_fasta(filename=seqb)

    # Generate dotplot
    print(colored(f'Generating DotPlot Matrix...', 'yellow'))
    dp = dotplot(seqA=seqA_string, seqB=seqB_string, w=int(w), s=int(s))

    # Generate results files
    if output.endswith(".txt"):
        # Print out in textfile as ASCII
        dotplot2Ascii(dp=dp,seqA=seqA_string, seqB=seqB_string, heading=title,
                      filename=output)

    if output.endswith(('.png','.ps','.pdf')):
        # Create and save graph file
        dotplot2Graphics(dp=dp, labelA=Path(seqa).name,labelB=Path(seqb).name,
                         heading=title,filename=output,spFlag=True, bigFlag=True)
    print(colored(f'Results saved: {output}', 'green'))

    return



if __name__ == '__main__':
    main()

# main(w=10, s=6, seqA='handout 03/human_pax6_NM_001604.fasta',
#      seqB='handout 03/mouse_pax6_NM_013627.fasta', title='human-vs-mouse',
#      output='human-vs-mouse.png')



