from Bio import SeqIO
from mimetypes import guess_type
from functools import partial
import parasail, time, types, gzip

def parseFasta(fi):
    sequences = []
    encoding = guess_type(fi)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(fi) as fasta:
        try:
            for record in SeqIO.parse(fasta, "fasta"):
                sequences.append(record)
        except:
            print('File format incorrect or file content does not correspond to fasta.')
    return sequences

def parseFastq(fi, suffix_array, k, ref):
    cmax_score = 0
    res = {}
    encoding = guess_type(fi)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(fi) as fastq:
        for record in SeqIO.parse(fastq, "fastq"):
            if str(record.seq[len(record.seq)//2:(len(record.seq)//2)+5]) != "NNNNN": 
                kmers = find_filtered_kmers(record.seq,k)
                seeds = find_Seeds(kmers, ref, suffix_array, k)
                for seed in seeds.values():
                    for elem in seed:
                        end_pos = min(elem+k+(len(record.seq)//2), len(ref))
                        start_pos = max(elem-(len(record.seq)//2), 0)
                        alignment = parasail.sg_dx_trace_scan(str(record.seq), str(ref[start_pos:end_pos]), 10, 1, parasail.dnafull)
                        #print(alignment.similar)
                        #print(alignment.traceback.query+"\n"+alignment.traceback.comp+"\n"+alignment.traceback.ref)
                        #alignment = parasail.sg_dx_trace_scan_sat(str(record.seq), str(ref), 10, 1, parasail.dnafull)
                        
                        if alignment.score >= 150 and alignment.score >= max_score:
                            max_score = alignment.score
                            max_alignment = alignment
                if max_score >= 150:
                    res[record.name] = [max_alignment, record.seq]
                    max_score = 0
    return res

def find_kmers(seq, k):
    kmers = []
    rc = reverseComplement(seq)
    for i in range(len(seq) - k+1):
        kmers.append(rc[i:i+k])
        kmers.append(seq[i:i+k])
    return kmers

def find_filtered_kmers(seq, k):
    kmers = []
    i = 0
    rc = reverseComplement(seq)
    while i <= len(seq) - k+1:
        kmers.append(rc[i:i+k])
        kmers.append(seq[i:i+k])
        i += 5
    return kmers

def reverseComplement(seq):
        pairs = {"A" : "T",
                 "T" : "A",
                 "C" : "G",
                 "G" : "C",
                 "N" : "N"}
        rc = ""
        for i in seq:
            rc += pairs[i]  #Pour chaque nucléotide, insérer son complément à la suite du complément des précédentes.
        return rc[::-1]

def create_suffix_array(sequence, k):
    suffix_array = list(range(len(sequence)-k))
    suffix_array.sort(key = lambda i: sequence[i:])
    return suffix_array

def dichotomicSearch(kmer, sequence, suffix_array, k):
    first = 0
    last = len(suffix_array)-1
    while first <= last:
        mid = first + (last - first) // 2
        if sequence[suffix_array[mid]:suffix_array[mid]+k] == kmer and sequence[suffix_array[mid-1]:suffix_array[mid-1]+k] < kmer:
            first = mid
            break
        elif sequence[suffix_array[mid]:suffix_array[mid]+k] < kmer:
            first = mid + 1
        else:
            last = mid -1

    start = first
    last = len(suffix_array)-1
    while start <= last:
        mid = start + (last - start) // 2
        if sequence[suffix_array[mid]:suffix_array[mid]+k] == kmer and mid == len(suffix_array) -1:
            last = mid
            break
        elif sequence[suffix_array[mid]:suffix_array[mid]+k] == kmer and sequence[suffix_array[mid+1]:suffix_array[mid+1]+k] > kmer:
            last = mid
            break
        elif sequence[suffix_array[mid]:suffix_array[mid]+k] < kmer:
            start = mid + 1
        else:
            last = mid -1
    return (first, last)

def find_Seeds(kmers, sequence, suffix_array, k):
    seeds = {}
    for j in range(len(kmers)):
        kmer = kmers[j]
        first, last = dichotomicSearch(kmer, sequence, suffix_array, k)
        for i in range(first, last+1):
            if(seeds.get(kmer, 0) != 0):    
                seeds[kmer].append(suffix_array[i])
            else:
                seeds[kmer] = [suffix_array[i]]
    
    return seeds

def write_output(filename, res, exec_time):
    with open(filename, 'w') as f:
        print(f"Total execution time: {exec_time}", file = f)
        for key in res.keys():
            cigar = res[key][0].cigar.decode
            name = key
            start_pos = res[key][2]
            score = res[key][0].score
            seq = res[key][1]
            print(name, start_pos, score, seq, cigar, sep='\t', file=f)
    