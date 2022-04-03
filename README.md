# Seed-Extend
Project consisting in comparing query sequences with reference genomes.
## Usage
seedExtend.py [-h] [-s] [-t] genome reads [reads ...] out k [gap]  

### positional arguments:
<font color="red">genome</font> Fasta file containing reference sequence  
<font color="red">reads</font>  Fasta file(s) containing query sequences  
<font color="red">out</font>    Output file  
<font color="red">k</font>     Word size of the kmers the alignment should use  
<font color="red">gap</font>   Gap between first nucleotides of consecutive kmers (1 by default)  
  
### optional arguments: 
<font color="red">-h</font>  show this help message and exit  
<font color="red">-s</font> Run the program on the first 1000 reads only  
<font color="red">-t</font>  Type of analysis to perform (S for seed-extend or F for filtering, S by default)  

### Examples:
`./SeedExtend.py Data/GCF_000851145.1_ViralMultiSegProj14892_genomic.fna Data/SRR10971381_1.fastq.gz result_grippe.txt 31`  
`./SeedExtend.py Data/GCF_000851145.1_ViralMultiSegProj14892_genomic.fna.gz Data/SRR10971381_1.fastq result_grippe.txt 15 15 -s`