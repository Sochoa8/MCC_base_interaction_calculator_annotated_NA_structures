# MCC_base_interaction_calculator_annotated_NA_structures
This Python script calculates Matthews Correlation Coefficient (MCC) and other performance metrics for nucleic acid structures. It processes RNApdbee output files to evaluate canonical, non-canonical, and stacking interactions, comparing predicted vs. reference structures.

The script uses the following output fies from the RNApdbee webserver (http://rnapdbee.cs.put.poznan.pl)

- Canonical base pairs (bpseq files)
- Non-canonical interactions (CSV files)
- Stacking interactions (CSV files)
  
##  Installation
This script requires **Python 3.1**. Clone the repository to get started:
```bash
git clone https://github.com/Sochoa8/MCC_base_interaction_calculator_annotated_NA_structures.git
cd MCC_base_interaction_calculator_annotated_NA_structures
```
## Usage
Run the script:
```
python MCC_scores.py
```
- Select (multiple) RNApdbee bpseq and CSV files when prompted.
- The script will extract base interactions and compute MCC scores.
- Results are saved in a structured CSV file.

## Expected Input Files
The script expects RNApdbee output files named using this pattern:

```filtered_<structure_id>_(af|pdb)[_-](2D-bpseq|non-canonical|stacking).<extension>```


### Example using PDB ID: 7D7W

For output files retrieved from MCannotate (hosted on RNApdbee) analysis of PDB structure

- Canonical base pairs (bpseq file):
  ``` filtered_7d7w_pdb-2D-bpseq.txt```
- Non-canonical interactions (CSV file):
  ```filtered_7d7w_pdb-non-canonical.csv```
- Stacking interactions (CSV files):
  ```filtered_7d7w_pdb-stacking.csv```

For output files retrieved from MCannotate (hosted on RNApdbee) analysis of AF predicted structure

- Canonical base pairs (bpseq file):
  ``` filtered_7d7w_af-2D-bpseq.txt```
- Non-canonical interactions (CSV file):
  ```filtered_7d7w_af-non-canonical.csv```
- Stacking interactions (CSV files):
  ```filtered_7d7w_af-stacking.csv```

### Input File Descriptions
Canonical Interactions (```2D-bpseq.txt```)

- Format: BPSEQ
- Example Structure:
```
1 G  0
2 C  12
3 A  0
4 U  0
5 G  10
6 C  0
7 G  0
8 G  0
9 G  0
10 C 5
11 A 0
12 G 2
```
- Column 1: Nucleotide index
- Column 2: Nucleotide type
- Column 3: Paired index (0 = unpaired)
  
Non-Canonical Interactions (```non-canonical.csv```)

- Format: Semicolon (;) delimited CSV
- Required Column: ```Base-pair``` (formatted as ```NucleotideX - NucleotideY```)
- Example:
```
Base-pair; Interaction Type
A.g1 - A.g6; Hoogsteen
G.c3 - U.c8; Sheared
```

Stacking Interactions (```stacking.csv```)

- Format: Semicolon (;) delimited CSV
- Required Column: ```Base-pair``` (``formatted as NucleotideX - NucleotideY``)
- Example:
```
Base-pair; Stack Strength
A.g1 - A.g2; Strong
G.c5 - C.c6; Moderate
```

##Example Output

The output is saved as a CSV file, summarizing MCC scores and other metrics.

| Structure | Residues | MCC (Canonical) | Precision (Canonical) | Recall (Canonical) | F1 (Canonical) | Accuracy (Canonical) | MCC (Non-Canonical) | Precision (Non-Canonical) | Recall (Non-Canonical) | F1 (Non-Canonical) | Accuracy (Non-Canonical) | MCC (Stacking) | Precision (Stacking) | Recall (Stacking) | F1 (Stacking) | Accuracy (Stacking) | MCC (Combined) | Precision (Combined) | Recall (Combined) | F1 (Combined) | Accuracy (Combined) |
|-----------|----------|----------------|----------------------|--------------------|--------------|--------------------|----------------------|------------------------|------------------------|------------------|----------------------|----------------|--------------------|------------------|--------------|------------------|--------------|--------------------|------------------|--------------|------------------|
| 148d      | 76       | 0.92           | 0.91                 | 0.94               | 0.92         | 0.95               | 0.85                 | 0.83                   | 0.88                   | 0.85             | 0.90                 | 0.78           | 0.80               | 0.76               | 0.78         | 0.84               | 0.88         | 0.88               | 0.89               | 0.88         | 0.92               |
| 1ehz      | 82       | 0.89           | 0.88                 | 0.90               | 0.89         | 0.93               | 0.81                 | 0.79                   | 0.85                   | 0.82             | 0.88                 | 0.75           | 0.77               | 0.74               | 0.76         | 0.82               | 0.86         | 0.87               | 0.87               | 0.86         | 0.90               |
| 2gdi      | 90       | 0.94           | 0.92                 | 0.95               | 0.93         | 0.96               | 0.87                 | 0.85                   | 0.90                   | 0.87             | 0.91                 | 0.80           | 0.82               | 0.78               | 0.80         | 0.86               | 0.90         | 0.91               | 0.91               | 0.90         | 0.93               |

## Notes:
- **MCC (Matthews Correlation Coefficient)**: Measures prediction quality. Values close to 1 indicate strong agreement between predicted and reference interactions.
- **Precision**: Fraction of correctly predicted interactions.
- **Recall**: Fraction of actual interactions correctly predicted.
- **F1-score**: Harmonic mean of precision and recall.
- **Accuracy**: Overall correctness of predictions.


  

