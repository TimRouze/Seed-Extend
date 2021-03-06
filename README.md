# Seed-Extend
Project consisting in comparing short reads from [SRR10971381](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR10971381) with reference genomes (Flu, Rhinovirus, HIV and Alpha Coronavirus). The goal is to associate the reads with a genome to understand from which organism the reads are from.
Part two of this project consists in comparing other short reads with [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) in order to filter the reads matching with the organism.
## Usage
SeedExtend.py [-h] [-g] [-s] [-t] genome out k reads [reads ...]
### positional arguments:
<font color="red">genome</font> Fasta file containing reference sequence  
<font color="red">reads</font>  Fasta file(s) containing query sequences  
<font color="red">out</font>    Output file  
<font color="red">k</font>     Word size of the kmers the alignment should use  
<font color="red">gap</font>   Gap between first nucleotides of consecutive kmers (1 by default)  
  
### optional arguments: 
<font color="red">-h</font>  show this help message and exit  
<font color="red">-s</font>  Stop the program once 1 seed meets the score threshold  
<font color="red">-t</font>  Type of analysis to perform (S for seed-extend or F for filtering, S by default)  

### Examples:
`./SeedExtend.py Data/GCF_000851145.1_ViralMultiSegProj14892_genomic.fna  result_grippe.txt 31 Data/SRR10971381_1.fastq.gz`  
`./SeedExtend.py Data/GCF_000851145.1_ViralMultiSegProj14892_genomic.fna.gz  result_grippe.txt 15 15 Data/SRR10971381_1.fastq -s`
