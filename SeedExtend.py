import sys, argparse
from Bio import SeqIO
from Utils import *
from mimetypes import guess_type
from functools import partial
import time, psutil, parasail, doctest

def seed_extend(files, suffix_array, k, ref, gap):
    """Core function which runs the seed-extend algorithm.

    Parameters
    ----------
    files : list of str
        A list of fasta files containing query sequences.
    suffix_array : list of int
        The suffix array of the reference sequence.
    k : int
        Word size of the kmers the alignment uses.
    ref : str
        Sequence of the reference sequence.
    gap : int
        _description_

    Returns
    -------
    res : dict of {str : tuple of (int, str, int)}
        _description_
    aux : tuple of (int, int)
        A tuple containing the memory used by the execution and the number of reads aligned.
    """
    cpt, cpt_extended, max_score, pos = 0, 0, 0, 0
    res = {}
    start = time.time()
    for fi in files:
        encoding = guess_type(fi)[1]
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
        with _open(fi) as fastq:
            for record in SeqIO.parse(fastq, "fastq"):
                if record.seq.count('N') <= len(record.seq)/2:
                    seeds = find_Seeds(ref, suffix_array, k, gap, record.seq)
                    for key in seeds.keys():
                        if key == "-":
                            record_seq = reverseComplement(record.seq)
                        else:
                            record_seq = record.seq
                        for elem in seeds[key]:
                            cpt_extended +=1
                            end_pos = min(elem+k+(len(record_seq)), len(ref))
                            start_pos = max(elem-(len(record_seq)), 0)
                            alignment = parasail.sg_dx_trace_scan(str(record_seq), str(ref[start_pos:end_pos]), 10, 1, parasail.dnafull)
                            #print(alignment.similar)
                            #print(alignment.traceback.query+"\n"+alignment.traceback.comp+"\n"+alignment.traceback.ref)
                            #alignment = parasail.sg_dx_trace_scan_sat(str(record.seq), str(ref), 10, 1, parasail.dnafull)
                            
                            if alignment.score >= 150 and alignment.score >= max_score:
                                max_score = alignment.score
                                pos = elem
                                max_alignment = alignment
                    if max_score >= 150:
                        res[record.name] = (max_alignment, record.seq, pos)
                        max_score, pos = 0, 0
                if cpt%1000 == 0:
                    curr = time.time() - start
                    sys.stdout.write(f"\r{((cpt/28282964)*100):.2f}% of total reads have been treated. computation time is {curr:.2f}secs")
                    sys.stdout.flush()
                cpt += 1
                    
    memory = psutil.Process().memory_info().rss / (1024 * 1024)
    aux = (memory, cpt_extended)
    return (res, aux)

def keep_matching_reads(files, suffix_array, k, ref, gap):
    cpt = 0
    res = []
    start = time.time()
    for fi in files:
        with gzip.open(fi, "rt") as fastq:
            for record in SeqIO.parse(fastq, "fastq"):
                if record.seq.count('N') <= len(record.seq)/2:
                    seeds = find_Seeds(ref, suffix_array, k, gap, record.seq)
                    for key in seeds.keys():
                        if key == "-":
                            record_seq = reverseComplement(record.seq)
                        else:
                            record_seq = record.seq
                        for elem in seeds[key]:
                            cpt_extended +=1
                            end_pos = min(elem+k+(len(record_seq)), len(ref))
                            start_pos = max(elem-(len(record_seq)), 0)
                            alignment = parasail.sg_dx_trace_scan(str(record_seq), str(ref[start_pos:end_pos]), 10, 1, parasail.dnafull)
                            #print(alignment.similar)
                            #print(alignment.traceback.query+"\n"+alignment.traceback.comp+"\n"+alignment.traceback.ref)
                            #alignment = parasail.sg_dx_trace_scan_sat(str(record.seq), str(ref), 10, 1, parasail.dnafull)
                            print(cpt)
                            if alignment.score >= 285:
                                res.append(record)
                                break
                        else:
                            continue
                        break
                if cpt%1000000 == 0:
                    curr = time.time() - start
                    sys.stdout.write(f"\r{((cpt/1000000)*100):.2f}% of total reads have been treated. computation time is {curr:.2f}secs")
                    sys.stdout.flush()
                    break
                cpt += 1
    return res


def main():
    #Command line options
    parser = argparse.ArgumentParser(description = 'Seed and extend alignment tool.')
    parser.add_argument('-g', '--genome', required=True, metavar='', help='Fasta file containing reference sequence')
    parser.add_argument('-r', '--reads', required=True, metavar='', nargs= '+', help='Fasta file(s) containing query sequences')
    parser.add_argument('-o', '--out', required=True, metavar='', help='Output file')
    parser.add_argument('-k', '--kmersize', required=True, metavar='', help='Word size of the kmers the alignment should use')
    parser.add_argument('-b', '--breaksize', required=False, metavar='', help='Gap between selected kmers (1 means no gap)', default=1)
    parser.add_argument('-t', '--type', required=False, metavar='', choices=['S', 's', 'F', 'f'], help='Type of analysis to perform (S for seed-extend or F for filtering)', default='S')
    args = parser.parse_args(sys.argv[1:])

    start = time.time()
    sequence = parseFasta(args.genome)
    suffix_array = create_suffix_array(sequence, int(args.kmersize))

    if args.type == 'S' or args.type =='s':
        #Seed-extend algorithm
        res, aux = seed_extend(args.reads, suffix_array, int(args.kmersize), sequence, int(args.breaksize))
        end = time.time()
        write_output(args.out, res, (end-start), aux)
    else:
        #filtering algorithm
        res = keep_matching_reads(args.reads, suffix_array, int(args.kmersize), sequence, int(args.breaksize))
        end = time.time()
        write_output(args.out, res, (end-start), [])
        

        
        print(f"\nExecution complete, please find the results in the {args.out} file.")
    """except Exception as e:
        with open("res_log.txt", "w") as err_file:
            print(e, file = err_file)"""
    

if __name__=="__main__":
    main()