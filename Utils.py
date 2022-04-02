from Bio import SeqIO
from mimetypes import guess_type
from functools import partial
import gzip

def parseFasta(fi):
    """Read a fasta or fasta.gz file and extract the sequence within.

    Parameters
    ----------
    fi : str
        Path of a fasta or fasta.gz file.

    Returns
    -------
    str
        The sequence in the input file.
    """
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

def reverseComplement(seq):
    """Compute the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        The sequence to reverse.

    Returns
    -------
    str
        The input sequence's reverse complement.

    Examples
    --------
    >>>reverseComplement('AGTTC')
    'GAACT'
    >>>reverseComplement('TCGACCTTCTCGTGGGAGGT')
    'ACCTCCCACGAGAAGGTCGA'
    """
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
    """Compute the suffix array of a DNA sequence.

    Parameters
    ----------
    sequence : str
        The DNA sequence we want the suffix array of.
    k : int
        The word size we are working with.

    Returns
    -------
    list of int
        The sorted suffix array of the input sequence, excluding suffixes of length < k.

    Examples
    --------
    >>>create_suffix_array('banana', 0)
    [6, 5, 3, 1, 0, 4, 2]
    >>>create_suffix_array('banana', 2)
    [3, 1, 0, 4, 2]
    >>>create_suffix_array('TCGACCTTCTCGTGGGAGGT', 3)
    [12, 0, 7, 10, 13, 6, 5, 4, 1, 17, 8, 2, 18 11, 9, 14, 15, 3, 16]
    >>>create_suffix_array('TCGACCTTCTCGTGGGAGGT', 12)
    [0, 7, 6, 5, 4, 1, 8, 2, 3]
    """
    suffix_array = list(range(len(sequence)-k+1))
    suffix_array.sort(key = lambda i: sequence[i:])
    return suffix_array

def dichotomicSearch(kmer, sequence, suffix_array, k):
    """Search for the suffixes in a sequences suffix array which start with a given kmer.

    Parameters
    ----------
    kmer : str
        A word of size k on the following alphabet: {A,C,T,G}
    sequence : str
        The DNA sequence we are working on.
    suffix_array : list of int
        The suffix array of the sequence.
    k : int
        The word size of kmer.

    Returns
    -------
    First: int
        The index of the first suffix in the suffix array which starts with the word kmer.
    Last: int
        The index of the last suffix in the suffix array which starts with the word kmer.

    Examples
    --------
    >>>dichtomicSearch('TT', 'AGTTC', [0, 1, 3, 2], 2)
    (4,4)
    >>>dichotomicSearch('ACG', 'ACGCTCCCACGACGAACGGTCGA', [14, 11, 8, 0, 15, 7, 6, 5, 20, 12, 9, 1, 16, 3, 13, 10, 2, 17, 18, 4, 19], 3)
    (2, 5)
    """
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
    """Compute all the seeds in a sequence which match with a kmer present in a given read.

    Parameters
    ----------
    seq : str
        The DNA sequence we are looking for the seeds in.
    suffix_array : list of int
        The suffix_array of seq.
    k : int
        The word size used to match kmers with seeds.
    gap : int
        Gap between the first nucleotide of each consecutive kmer to be selected from tne read.
    read : str
        The DNA sequence we take kmers from.

    Returns
    -------
    dict of {str : int}
        A dicitonary which maps suffix indexes which start with a kmer present in the read with the strand (- or +) that seed is on.
    """
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
                seeds["+"].append(suffix_array[j])
            else:
                seeds["+"] = [suffix_array[j]]
        for j in range(first_rc, last_rc+1):
            if(seeds.get(kmer_rc, 0) != 0):    
                seeds["-"].append(suffix_array[j])
            else:
                seeds["-"] = [suffix_array[j]]
    return seeds

def write_output(filename, res, exec_time, aux):
    """Write the results of the execution of the program.

    Parameters
    ----------
    filename : str
        Path of the output file results should be written in.
    res : dict of {str : tuple of (int, str, int)}
        _description_
    exec_time : int
        Time the execution took.
    aux : tuple of (int, int)
        A tuple containing the memory used by the execution and the number of reads aligned.
    """
    with open(filename, 'w') as f:
        if aux == []:
            print(f"Total execution time: {exec_time}", file = f)
            for elem in res:
                print(elem.id, elem.name, elem.seq, sep='\t', file=f)
        else:
            print(f"Total execution time: {exec_time}, {aux[1]} reads where aligned and {aux[0]}Mb of memory was used", file = f)
            for key in res.keys():
                cigar = res[key][0].cigar.decode
                name = key
                start_pos = res[key][2]
                score = res[key][0].score
                seq = res[key][1]
                print(name, start_pos, score, seq, cigar, sep='\t', file=f)
    
### UNUSED FUNCTIONS ###

def find_kmers(seq, k):
    """Read a DNA sequence and extract all the words of size k in it.

    Parameters
    ----------
    seq : str
        A DNA sequence.
    k : int
        The word size used.

    Returns
    -------
    list of str
        A list of all the words of size k found in the sequence.
    """
    kmers = []
    for i in range(len(seq) - k+1):
        kmers.append(seq[i:i+k])
    return kmers

def find_filtered_kmers(seq, k):
    """Read a DNA sequence and extract the consecutive words of size k in it.

    Parameters
    ----------
    seq : str
        A DNA sequence.
    k : int
        The word size used.

    Returns
    -------
    list of str
        A list of the consecutive words size k found in the sequence.
    """
    kmers = []
    i = 0
    while i <= len(seq) - k+1:
        kmers.append(seq[i:i+k])
        i += k
    return kmers