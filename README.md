# Smith-waterman-alignment

This is a python script for a smith waterman algorithm. 
It aligns two sequences by matches/mismatches (also known as substitutions), insertions, and deletions.
Both insertions and deletions are the operations that introduce gaps, which are represented by dashes.
The code creates a matrix of the scores for each match/mismatch and then creates a traceback matrix to find the most parsimonious sequence alignment.
Pentalties and match/mismatch scoreds can be altered as they are stored as variables.
