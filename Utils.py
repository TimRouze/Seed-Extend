from Bio import SeqIO
from mimetypes import guess_type
from functools import partial
import gzip

def parseFasta(fi):
    sequence = ""
    encoding = guess_type(fi)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(fi) as fasta:
        try:
            for record in SeqIO.parse(fasta, "fasta"):
                sequence += record.seq
        except:
            print('File format incorrect or file content does not correspond to fasta.')
    return sequence

def find_kmers(seq, k):
    kmers = []
    for i in range(len(seq) - k+1):
        kmers.append(seq[i:i+k])
    return kmers

def find_filtered_kmers(seq, k):
    kmers = []
    i = 0
    #rc = reverseComplement(seq)
    while i <= len(seq) - k+1:
        #kmers.append(rc[i:i+k])
        kmers.append(seq[i:i+k])
        i += k
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

def find_Seeds(seq, suffix_array, k, gap, read):
    seeds = {}
    kmers, kmers_rc = [], []
    i = 0
    rc = reverseComplement(read)
    while i <= len(read) - k+1:
        kmer_rc = rc[i:i+k]
        kmer = read[i:i+k]
        i += gap
        first, last = dichotomicSearch(kmer, seq, suffix_array, k)
        first_rc, last_rc = dichotomicSearch(kmer_rc, seq, suffix_array, k)
        for j in range(first, last+1):
            if(seeds.get(kmer, 0) != 0):
                print(first, last)
                seeds["normal"].append(suffix_array[j])
            else:
                seeds["normal"] = [suffix_array[j]]
        for j in range(first_rc, last_rc+1):
            if(seeds.get(kmer_rc, 0) != 0):    
                seeds["rc"].append(suffix_array[j])
            else:
                seeds["rc"] = [suffix_array[j]]
    return seeds

def write_output(filename, res, exec_time, aux):
    with open(filename, 'w') as f:
        if len(aux) == len(res):
            print(f"Total execution time: {exec_time}", file = f)
            for i in range(len(res)):
                print(res[i].id, res[i].seq, aux[i], sep='\t', file=f)
        else:
            print(f"Total execution time: {exec_time}, {aux[1]} reads where aligned and {aux[0]}Mb of memory was used", file = f)
            for key in res.keys():
                cigar = res[key][0].cigar.decode
                name = key
                start_pos = res[key][2]
                score = res[key][0].score
                seq = res[key][1]
                print(name, start_pos, score, seq, cigar, sep='\t', file=f)
    