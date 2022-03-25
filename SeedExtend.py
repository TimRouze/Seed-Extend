import sys, argparse
from Bio import SeqIO
from Utils import *
import parasail, time

def find_Seeds(kmers, sequence, suffix_array, k):
    seeds = {}
    for j in range(len(kmers)):
        kmer = kmers[j]
        first = dichotomicSearchFirst(kmer, sequence, suffix_array, k)
        last = dichotomicSearchLast(kmer, sequence, suffix_array, k)
        for i in range(first, last+1):
            if(seeds.get(kmer, 0) != 0):    
                seeds[kmer].append(suffix_array[i])
            else:
                seeds[kmer] = [suffix_array[i]]
    return seeds


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
    queries = parseFastq(args.reads)
    suffix_array = create_suffix_array(sequences[0].seq, int(args.kmersize))
    for query in queries:
        kmers = find_kmers(query.seq,int(args.kmersize))
        seeds = find_Seeds(kmers, sequences[0].seq, suffix_array, int(args.kmersize))
        if (len(seeds) > 0):
            alignment = parasail.sg_dx_trace_scan(str(query.seq), str(sequences[0].seq), 10, 1, parasail.dnafull)
            #print(alignment.score)
            #print(alignment.traceback.query+"\n"+alignment.traceback.comp+"\n"+alignment.traceback.ref)
            alignment.cigar.decode
            #alignment = parasail.sg_dx_trace_scan_sat(str(query.seq), str(sequences[0].seq), 10, 1, parasail.dnafull)
            #print(seeds)
    end = time.time()
    print(end - start)

if __name__=="__main__":
    main()