# Seed-Extend
Project consisting in comparing query sequences with reference genomes.
## Usage
SeedExtend.py [-h] -g  -r  [...] -o  -k  [-b] [-t]

optional arguments:  
-h, show this help message and exit  
-g, Fasta file containing reference sequence  
-r, Fasta file(s) containing query sequences  
-o, Output file  
-k, Word size of the kmers the alignment should use  
-b, Gap between selected kmers (1 means no gap)  
-t, Type of analysis to perform (S for seed-extend or F for filtering)  
