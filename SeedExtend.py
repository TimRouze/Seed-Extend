import sys, argparse
from Bio import SeqIO
from Utils import *
import time, psutil



def main():
    #Command line options
    parser = argparse.ArgumentParser(description = 'Simple ORF finder.')
    parser.add_argument('-g', '--genome', required=True, metavar='', help='Fasta containing reference sequence')
    parser.add_argument('-r', '--reads', required=True, metavar='', nargs= '+', help='Fasta containing query sequences')
    parser.add_argument('-o', '--out', required=True, metavar='', help='Output file')
    parser.add_argument('-k', '--kmersize', required=True, metavar='', help='size of the kmers')
    args = parser.parse_args(sys.argv[1:])
    start = time.time()
    sequence = parseFasta(args.genome)
    suffix_array = create_suffix_array(sequence, int(args.kmersize))
    try:
        res, aux = parseFastq(args.reads, suffix_array, int(args.kmersize), sequence)
        end = time.time()

        write_output(args.out, res, (end-start), aux)
        print(f"\nExecution complete, please find the results in the {args.out} file.")
    except Exception as e:
        with open("res_log.txt", "w") as err_file:
            print(e, file = err_file)
            if hasattr(e, 'message'):
                print(e.message, file = err_file)
            elif hasattr(e, 'strerror'):
                print(e.strerror, file = err_file)
            else:
                print(e, file = err_file)
    

if __name__=="__main__":
    main()