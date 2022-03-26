from Bio import SeqIO
import parasail

def parseFasta(fi):
    sequences = []
    with open(fi, "r") as fasta:
        try:
            for record in SeqIO.parse(fasta, "fasta"):
                sequences.append(record)
        except:
            print('File format incorrect or file content does not correspond to fasta.')
    return sequences

def parseFastq(fi):
    sequences = []
    with open(fi, "r") as fastq:
        try:
            for record in SeqIO.parse(fastq, "fastq"):
                sequences.append(record)
                if len(sequences) == 1000:
                    break;
        except:
            print("File format incorrect or file content does not correspond to fastq.")
    return sequences

def parseFastq2(fi, suffix_array, k, ref):
    sequences = []
    with open(fi, "r") as fastq:
        try:
            for record in SeqIO.parse(fastq, "fastq"):
                kmers = find_kmers(record.seq,k)
                seeds = find_Seeds(kmers, ref, suffix_array, k)
                if (len(seeds) > 0):
                    alignment = parasail.sg_dx_trace_scan(str(record.seq), str(ref), 10, 1, parasail.dnafull)
                    #print(alignment.score)
                    #print(alignment.traceback.query+"\n"+alignment.traceback.comp+"\n"+alignment.traceback.ref)
                    alignment.cigar.decode
                    #alignment = parasail.sg_dx_trace_scan_sat(str(query.seq), str(sequences[0].seq), 10, 1, parasail.dnafull)
                    #print(seeds)
                sequences.append(record)
                if len(sequences) == 1000:
                    break;
        except:
            print("File format incorrect or file content does not correspond to fastq.")
    return sequences

def find_kmers(seq, k):
    kmers = []
    i = 0
    while i < len(seq) - k:
        if len(seq[i:i+k]) == k: #and sequence[i:i+k] not in kmers:
            kmers.append(seq[i:i+k]) 
        i += k
    return kmers

def reverseComplement(seq):
        pairs = {"A" : "T",
                 "T" : "A",
                 "C" : "G",
                 "G" : "C"}
        rc = ""
        for i in seq:
            rc += pairs[i]  #Pour chaque nucléotide, insérer son complément à la suite du complément des précédentes.
        return rc[::-1]

def create_suffix_array(sequence, k):
    #sequence.append("$")
    suffix_array = list(range(len(sequence)-k))
    suffix_array.sort(key = lambda i: sequence[i:])
    return suffix_array

def dichotomicSearchFirst(kmer, sequence, suffix_array, k):
    first = 0
    last = len(suffix_array)-1
    while first <= last:
        mid = first + (last - first) // 2
        if sequence[suffix_array[mid]:suffix_array[mid]+k] == kmer and sequence[suffix_array[mid-1]:suffix_array[mid-1]+k] < kmer:
            return mid
        elif sequence[suffix_array[mid]:suffix_array[mid]+k] < kmer:
            first = mid + 1
        else:
            last = mid - 1
    return first

def dichotomicSearchLast(kmer, sequence, suffix_array, k):
    first = 0
    last = len(suffix_array)-1
    while first <= last:
        mid = first + (last - first) // 2
        if sequence[suffix_array[mid]:suffix_array[mid]+k] == kmer and sequence[suffix_array[mid+1]:suffix_array[mid+1]+k] > kmer:
            return mid
        elif sequence[suffix_array[mid]:suffix_array[mid]+k] < kmer:
            first = mid + 1
        else:
            last = mid - 1
    return last


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
