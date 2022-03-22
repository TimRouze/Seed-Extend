import sys, argparse
from Bio import SeqIO
from Utils import *

def find_Seeds(kmers, sequence, suffix_array, k):
    seeds = {}
    for kmer in kmers:
        first = dichotomicSearchFirst(kmer, sequence, suffix_array, k)
        last = dichotomicSearchLast(kmer, sequence, suffix_array, k)
        for i in range(first, last+1):
            seeds[kmer].append(suffix_array[i])
    return seeds

def main():
    #Command line options
    parser = argparse.ArgumentParser(description = 'Simple ORF finder.')
    parser.add_argument('-g', '--genome', required=True, metavar='', help='Fasta containing reference sequence')
    parser.add_argument('-r', '--reads', required=True, metavar='', help='Fasta containing query sequences')
    parser.add_argument('-o', '--out', required=True, metavar='', help='Output file')
    parser.add_argument('-k', '--kmersize', required=True, metavar='', help='size of the kmers')
    args = parser.parse_args(sys.argv[1:])
    sequences = parseFasta(args.genome)
    queries = parseFasta(args.reads)
    kmers = find_kmers(queries[0],args.kmersize)
    suffix_array = create_suffix_array(sequences[0].seq, args.kmersize)
    print(len(suffix_array))
    print(len(sequences[0].seq))

if __name__=="__main__":
    main()