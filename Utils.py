from ctypes import alignment
from Bio import SeqIO
import parasail, time

def parseFasta(fi):
    sequences = []
    with open(fi, "r") as fasta:
        try:
            for record in SeqIO.parse(fasta, "fasta"):
                sequences.append(record)
        except:
            print('File format incorrect or file content does not correspond to fasta.')
    return sequences

def parseFastq(fi, suffix_array, k, ref):
    cpt, seed = 0, 0
    records, alignements = [], []
    with open(fi, "r") as fastq:
        try:
            for record in SeqIO.parse(fastq, "fastq"):
                if str(record.seq[len(record.seq)//2:(len(record.seq)//2)+5]) != "NNNNN":
                    kmers = find_filtered_kmers(record.seq,k)
                    seeds = find_Seeds(kmers, ref, suffix_array, k)
                    if (len(seeds) > 0):
                        #cpt += 1
                        alignment = parasail.sg_dx_trace_scan(str(record.seq), str(ref), 10, 1, parasail.dnafull)
                        #print(alignment.similar)
                        #print(alignment.traceback.query+"\n"+alignment.traceback.comp+"\n"+alignment.traceback.ref)
                        #alignment = parasail.sg_dx_trace_scan_sat(str(record.seq), str(ref), 10, 1, parasail.dnafull)
                        score = alignment.score
                        if score >= 150:
                            records.append(record)
                            alignements.append(alignment)
                    #if cpt == 1000:
                        #break
        except:
            print("File format incorrect or an error occured during seeding.")
    return (records, alignements)

def find_kmers(seq, k):
    kmers = []
    for i in range(len(seq) - k+1):
        kmers.append(seq[i:i+k])
    return kmers

def find_filtered_kmers(seq, k):
    kmers = []
    i = 0
    while i <= len(seq) - k+1:
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
        """if sequence[suffix_array[mid]:suffix_array[mid]+k] == kmer and mid == len(suffix_array) -1:
            return mid"""
        if sequence[suffix_array[mid]:suffix_array[mid]+k] == kmer and sequence[suffix_array[mid+1]:suffix_array[mid+1]+k] > kmer:
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

def write_output(filename, recs, aligns):
    with open(filename, 'w') as f:
        for i in range(len(recs)):
            cigar = aligns[i].cigar.decode
            name = recs[i].name
            for i in range(len(cigar)):
                if cigar[i] == "D":
                    start_pos = cigar[0:i]
                    break
            score = aligns[i].score
            seq = recs[i].seq
            print(name, start_pos, score, seq, cigar, sep='\t', file=f)
    