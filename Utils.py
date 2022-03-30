from Bio import SeqIO
import gzip
import parasail, time, sys, psutil

def parseFasta(fi):
    sequence = ""
    with gzip.open(fi, "rt") as fasta:
        try:
            for record in SeqIO.parse(fasta, "fasta"):
                sequence += record.seq
        except:
            print('File format incorrect or file content does not correspond to fasta.')
    return sequence

def parseFastq(files, suffix_array, k, ref):
    cpt, cpt_extended, max_score, pos = 0, 0, 0, 0
    res = {}
    start = time.time()
    for fi in files:
        with gzip.open(fi, "rt") as fastq:
            for record in SeqIO.parse(fastq, "fastq"):
                if record.seq.count('N') <= len(record.seq)/2:
                    kmers = find_filtered_kmers(record.seq,k)
                    seeds = find_Seeds(kmers, ref, suffix_array, k)
                    for seed in seeds.values():
                        cpt_extended += 1
                        for elem in seed:
                            end_pos = min(elem+k+(len(record.seq)), len(ref))
                            start_pos = max(elem-(len(record.seq)), 0)
                            alignment = parasail.sg_dx_trace_scan(str(record.seq), str(ref[start_pos:end_pos]), 10, 1, parasail.dnafull)
                            #print(alignment.similar)
                            #print(alignment.traceback.query+"\n"+alignment.traceback.comp+"\n"+alignment.traceback.ref)
                            #alignment = parasail.sg_dx_trace_scan_sat(str(record.seq), str(ref), 10, 1, parasail.dnafull)
                            
                            if alignment.score >= 150 and alignment.score >= max_score:
                                max_score = alignment.score
                                pos = elem
                                max_alignment = alignment
                    if max_score >= 150:
                        res[record.name] = [max_alignment, record.seq, pos]
                        max_score, pos = 0, 0
                if cpt%1000 == 0:
                    curr = time.time() - start
                    sys.stdout.write(f"\r{((cpt/28282964)*100):.2f}% of total reads have been treated. computation time is {curr:.2f}secs")
                    sys.stdout.flush()
                cpt += 1
                    
    memory = psutil.Process().memory_info().rss / (1024 * 1024)
    aux = [memory, cpt_extended]
    return (res, aux)

def find_kmers(seq, k):
    kmers = []
    for i in range(len(seq) - k+1):
        kmers.append(seq[i:i+k])
    return kmers

def find_filtered_kmers(seq, k):
    kmers = []
    i = 0
    rc = reverseComplement(seq)
    while i <= len(seq) - k+1:
        kmers.append(rc[i:i+k])
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

def write_output(filename, res, exec_time, aux):
    with open(filename, 'w') as f:
        print(f"Total execution time: {exec_time}, {aux[1]} reads where aligned and {aux[0]}Mb of memory was used", file = f)
        for key in res.keys():
            cigar = res[key][0].cigar.decode
            name = key
            start_pos = res[key][2]
            score = res[key][0].score
            seq = res[key][1]
            print(name, start_pos, score, seq, cigar, sep='\t', file=f)
    