import sys, argparse
from Bio import SeqIO
from Utils import *
import time



def main():
    #Command line options
    parser = argparse.ArgumentParser(description = 'Simple ORF finder.')
    parser.add_argument('-g', '--genome', required=True, metavar='', help='Fasta containing reference sequence')
    parser.add_argument('-r', '--reads', required=True, metavar='', help='Fasta containing query sequences')
    parser.add_argument('-o', '--out', required=True, metavar='', help='Output file')
    parser.add_argument('-k', '--kmersize', required=True, metavar='', help='size of the kmers')
    args = parser.parse_args(sys.argv[1:])
    start = time.time()
    sequences = parseFasta(args.genome)
    suffix_array = create_suffix_array(sequences[0].seq, int(args.kmersize))
    recs, aligns = parseFastq(args.reads, suffix_array, int(args.kmersize), sequences[0].seq)
    end = time.time()
    print(end - start)
    write_output(args.out, recs, aligns)
    

if __name__=="__main__":
    main()