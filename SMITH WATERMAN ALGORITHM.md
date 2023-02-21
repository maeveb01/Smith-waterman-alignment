```python
#import packages
import numpy as np

#define sequences
seq1 = "CAGCCCCGTGGTATT"
seq2 = "CAGCGCCCGTGGTCT"
```


```python
#CREATE SCORING AND TRACEBACK MATRIX

def swa(seq1, seq2, MATCH =2, MISMATCH = -1, GAP_PENALTY = -2, STOP = 0, LEFT = 1, UP = 2, DIAGONAL =3):
    
    # MATRIX PARAMETERS - must start with row of zeros hence +1 when calculating length of matrix and -1 when finding cell with maximum value
    scoring_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1), np.int)  
    traceback_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1), np.int)  
    maximum_score = -1
    max_cell_index = (-1, -1)
    
    # CALCULATING HIGHEST SCORE FOR EVERY CELL IN MATRIX
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            
            # MATCH
            match_score = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            diagonal_score = scoring_matrix[i - 1, j - 1] + match_score
            # INDEL - verticle and horizonal mismated and need insertion in either row
            deletion = scoring_matrix[i - 1, j] + GAP_PENALTY
            insertion = scoring_matrix[i, j - 1] + GAP_PENALTY
            # HIGHEST SCORE
            scoring_matrix[i, j] = max(diagonal_score, deletion, insertion)
            
            # FINDING THE CELL WITH THE HIGHEST SCORE   
            if scoring_matrix[i, j] == 0: 
                traceback_matrix[i, j] = STOP     
            elif scoring_matrix[i, j] == insertion: 
                traceback_matrix[i, j] = LEFT      
            elif scoring_matrix[i, j] == deletion: 
                traceback_matrix[i, j] = UP      
            elif scoring_matrix[i, j] == diagonal_score: 
                traceback_matrix[i, j] = DIAGONAL     
            if scoring_matrix[i, j] >= maximum_score:
                max_cell_index = (i,j)
                maximum_score = scoring_matrix[i, j]
    
    # TRACEBACK MATRIX
    
    # TRACEBACK VARIABLES
    seq1_alignment = ""
    seq2_alignment = ""   
    gap_seq1 = ""   
    gap_seq2 = ""  
    (max_i, max_j) = max_cell_index
    
    #GOING BACK THROUGH THE MATRIX TO FIND NEXT HIGHEST SCORE
    while traceback_matrix[max_i, max_j] != STOP:
        if traceback_matrix[max_i, max_j] == DIAGONAL:
            gap_seq1 = seq1[max_i - 1]
            gap_seq2 = seq2[max_j - 1]
            max_i = max_i - 1
            max_j = max_j - 1
            
        elif traceback_matrix[max_i, max_j] == UP:
            gap_seq1 = seq1[max_i - 1]
            gap_seq2 = '-'
            max_i = max_i - 1    
            
        elif traceback_matrix[max_i, max_j] == LEFT:
            gap_seq1 = '-'
            gap_seq2 = seq2[max_j - 1]
            max_j = max_j - 1
            
        seq1_alignment = seq1_alignment + gap_seq1
        seq2_alignment = seq2_alignment + gap_seq2
    
    #CALCULATE NUMBER OF MATCHES, MISMATCHES, GAPS
    match = 0
    mismatch = 0 
    gaps = (seq1_alignment.count('-'))+(seq2_alignment.count('-'))

    for i in range(len(seq1_alignment)): 
        if seq1_alignment[i] == seq2_alignment[i]: 
            match += 1
        else: 
            mismatch += 1
    
    #PRINT MATRIX, ALIGNED SEQUENCES, AND STATS - reverse the sequences as they were organised from last column to first
    return print('Scoring Matrix:', '\n', str(scoring_matrix).replace('  [', '').replace('[', '').replace(']', ''), '\n', '\n', 'Traceback Matrix:', '\n',  str(traceback_matrix).replace('  [', ' ').replace('[', '').replace(']', ''), '\n', '\n',  'Aligned sequence:', '\n', seq1_alignment[::-1], '\n', seq2_alignment[::-1], '\n','\n', 'Score:', maximum_score, '\n', 'Matches:', match, '\n', 'Mismatches:', (mismatch-gaps), '\n', 'Gaps:', gaps)
```


```python
swa(seq1, seq2)
```

    Scoring Matrix: 
      0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
      0  2  0 -1  2  0  2  2  2  0 -1 -1 -1 -1  2  0
      0  0  4  2  0  1  0  1  1  1 -1 -2 -2 -2  0  1
      0 -1  2  6  4  2  0 -1  0  3  1  1  0 -2 -2 -1
      0  2  0  4  8  6  4  2  1  1  2  0  0 -1  0 -2
      0  2  1  2  6  7  8  6  4  2  0  1 -1 -1  1 -1
      0  2  1  0  4  5  9 10  8  6  4  2  0 -2  1  0
      0  2  1  0  2  3  7 11 12 10  8  6  4  2  0  0
      0  0  1  3  1  4  5  9 10 14 12 10  8  6  4  2
      0 -1 -1  1  2  2  3  7  8 12 16 14 12 10  8  6
      0 -1 -2  1  0  4  2  5  6 10 14 18 16 14 12 10
      0 -1 -2  0  0  2  3  3  4  8 12 16 20 18 16 14
      0 -1 -2 -2 -1  0  1  2  2  6 10 14 18 22 20 18
      0 -1  1 -1 -3 -2 -1  0  1  4  8 12 16 20 21 19
      0 -1 -1  0 -2 -4 -3 -2 -1  2  6 10 14 18 19 23
      0 -1 -2 -2 -1 -3 -5 -4 -3  0  4  8 12 16 17 21 
     
     Traceback Matrix: 
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
     0 3 0 3 3 0 3 3 3 0 3 3 3 3 3 0
     0 0 3 1 0 3 0 3 3 3 1 3 3 3 0 3
     0 3 2 3 1 1 0 2 0 3 1 3 0 1 2 2
     0 3 0 2 3 1 1 1 3 2 3 0 0 3 0 1
     0 3 3 2 2 3 3 1 1 1 0 3 1 3 3 1
     0 3 3 0 2 2 3 3 1 1 1 1 0 1 3 0
     0 3 3 0 2 2 2 3 3 1 1 1 1 1 0 0
     0 0 3 3 1 3 2 2 2 3 1 1 1 1 1 1
     0 3 2 2 3 2 2 2 2 2 3 1 1 1 1 1
     0 3 3 3 0 3 1 2 2 2 2 3 1 1 1 1
     0 3 3 0 0 2 3 2 2 2 2 2 3 1 1 1
     0 3 3 2 3 0 2 3 2 2 2 2 2 3 1 1
     0 3 3 1 1 2 2 0 3 2 2 2 2 2 3 1
     0 3 2 0 1 1 2 2 2 2 2 2 2 2 2 3
     0 3 3 2 3 1 1 2 2 0 2 2 2 2 2 2 
     
     Aligned sequence: 
     CAGC-CCCGTGGTAT 
     CAGCGCCCGTGGTCT 
     
     Score: 23 
     Matches: 13 
     Mismatches: 1 
     Gaps: 1

