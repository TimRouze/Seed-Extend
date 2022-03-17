from Bio import SeqIO

def parseFasta(fi):
    sequences = []
    with open(fi, "r") as fasta:
        try:
            for record in SeqIO.parse(fasta, "fasta"):
                sequences.append(record)
        except:
            print('File format incorrect or file content does not correspond to fasta.')
    return sequences

def kmers(seq, k):
    sequence = seq.seq
    kmers = []
    i = 0
    while i < len(sequence) - k:
        if len(sequence[i:i+k]) == k:
            kmers.append(str(sequence[i:i+k]))
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