```python
#import packages and define sequences
import numpy as np

seq1 = "CAGCCCCGTGGTATT"
seq2 = "CAGCGCCCGTGGTCT"
```


```python
#CREATE SCORING AND TRACEBACK MATRIX

MATCH = 2
MISMATCH = -1
GAP_PENALTY = -2

STOP = 0
LEFT = 1 
UP = 2
DIAGONAL = 3 
    
def scoring_matrix(seq1, seq2):
    
    # MATRIX PARAMETERS - must start with row of zeros hense +1
    scoring_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1), np.int)  
    traceback_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1), np.int)  
    
    # INITIALISE VARIABLES - ACCOUNT FOR ROW OF ZEROS
    max_score = -1
    max_index = (-1, -1)
    
    # CALCULATING SCORE FOR EVERY CELL IN MATRIX
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            
            # MATCH
            match_score = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            diagonal_score = scoring_matrix[i - 1, j - 1] + match_score
            
            # INDEL
            vertical_score = scoring_matrix[i - 1, j] + GAP_PENALTY
            horizontal_score = scoring_matrix[i, j - 1] + GAP_PENALTY
            
            # HIGHEST SCORE
            scoring_matrix[i, j] = max(diagonal_score, vertical_score, horizontal_score)
            
            # WHICH CELL HAS HIGHEST SCORE   
            if scoring_matrix[i, j] == 0: 
                traceback_matrix[i, j] = STOP
                
            elif scoring_matrix[i, j] == horizontal_score: 
                traceback_matrix[i, j] = LEFT
                
            elif scoring_matrix[i, j] == vertical_score: 
                traceback_matrix[i, j] = UP
                
            elif scoring_matrix[i, j] == diagonal_score: 
                traceback_matrix[i, j] = DIAGONAL 
                
            if scoring_matrix[i, j] >= max_score:
                max_index = (i,j)
                max_score = scoring_matrix[i, j]
    
    return scoring_matrix
```


```python
scoring_matrix(seq1, seq2)
```




    array([[ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
           [ 0,  2,  0, -1,  2,  0,  2,  2,  2,  0, -1, -1, -1, -1,  2,  0],
           [ 0,  0,  4,  2,  0,  1,  0,  1,  1,  1, -1, -2, -2, -2,  0,  1],
           [ 0, -1,  2,  6,  4,  2,  0, -1,  0,  3,  1,  1,  0, -2, -2, -1],
           [ 0,  2,  0,  4,  8,  6,  4,  2,  1,  1,  2,  0,  0, -1,  0, -2],
           [ 0,  2,  1,  2,  6,  7,  8,  6,  4,  2,  0,  1, -1, -1,  1, -1],
           [ 0,  2,  1,  0,  4,  5,  9, 10,  8,  6,  4,  2,  0, -2,  1,  0],
           [ 0,  2,  1,  0,  2,  3,  7, 11, 12, 10,  8,  6,  4,  2,  0,  0],
           [ 0,  0,  1,  3,  1,  4,  5,  9, 10, 14, 12, 10,  8,  6,  4,  2],
           [ 0, -1, -1,  1,  2,  2,  3,  7,  8, 12, 16, 14, 12, 10,  8,  6],
           [ 0, -1, -2,  1,  0,  4,  2,  5,  6, 10, 14, 18, 16, 14, 12, 10],
           [ 0, -1, -2,  0,  0,  2,  3,  3,  4,  8, 12, 16, 20, 18, 16, 14],
           [ 0, -1, -2, -2, -1,  0,  1,  2,  2,  6, 10, 14, 18, 22, 20, 18],
           [ 0, -1,  1, -1, -3, -2, -1,  0,  1,  4,  8, 12, 16, 20, 21, 19],
           [ 0, -1, -1,  0, -2, -4, -3, -2, -1,  2,  6, 10, 14, 18, 19, 23],
           [ 0, -1, -2, -2, -1, -3, -5, -4, -3,  0,  4,  8, 12, 16, 17, 21]])




```python
def traceback_matrix(H, b, b_='', old_i=0):
    # flip H to get index of **last** occurrence of H.max() with np.argmax()
    H_flip = np.flip(np.flip(H, 0), 1)
    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()
    if H[i, j] == 0:
        return b_, j
    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_
    return traceback(H[0:i, 0:j], b, b_, i)
```
