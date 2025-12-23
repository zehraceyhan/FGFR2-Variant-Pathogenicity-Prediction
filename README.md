# FGFR2-Variant-Pathogenicity-Prediction
Classification of FGFR2 protein mutations using protein sequences

# Overview
This project investigates the pathogenicity of missense variants in the FGFR2 gene by generating full-length mutant protein sequences and analyzing sequence-based features.

The reference protein sequence is obtained from Uniprot and known classifications of mutations are derived from Clinvar. 
When unknown variants are detected, the predictor will use BLOSUM62 matrix to calculate the mutation's value and calssify it based on its score. 

To ease the use of the predictor, mutation variants can be generated through the "mutantvar.py" file by randomly assigning an amino acid to a preffered position in the sequence. 

# How to run the predictor:
1. (not compulsory) generate a mutant variant sequence in "mutantvar.py" file 
2. place variant sequence as a string as input in line 216 of "predictor.py" file
3. run the predictor


