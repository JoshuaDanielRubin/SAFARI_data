import gzip
from typing import List, Dict, Tuple

# 1. Parsing the Fasta File
def parse_fasta(file_path: str) -> str:
    with open(file_path, 'r') as file:
        sequence = ''.join(line.strip() for line in file if not line.startswith('>'))
    return sequence

# 2. Compute Minimizer Index
def compute_minimizer(sequence: str, k: int, w: int) -> List[str]:
    if k > w:
        raise ValueError("k should be less than or equal to w")
    
    minimizers = []
    for i in range(len(sequence) - w + 1):
        window = sequence[i:i+w]
        minimizers.append(min([window[j:j+k] for j in range(w - k + 1)]))
    return minimizers

# 3. Rymer Transformation
def rymer_transform(sequence: str) -> str:
    return sequence.replace("T", "C").replace("G", "A")

# 4. Create Index Table
def create_index_table(sequence: str, k: int, w: int) -> Dict[str, List[int]]:
    index_table = {}
    minimizers = compute_minimizer(sequence, k, w)
    for i, minimizer in enumerate(minimizers):
        if minimizer not in index_table:
            index_table[minimizer] = []
        index_table[minimizer].append(i)
    return index_table

# 5. Reverse Complement
def reverse_complement(seq: str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq))

# 6. Parse Gzipped Fastq
def parse_gzipped_fastq(file_path: str) -> List[str]:
    sequences = []
    with gzip.open(file_path, 'rt') as file:
        lines = file.readlines()
        for i in range(1, len(lines), 4):
            sequences.append(lines[i].strip())
    return sequences

# 7. Count Matches
def count_matches(reads: List[str], k: int, w: int, minimizer_table: Dict[str, List[int]], rymer_table: Dict[str, List[int]]) -> Dict[str, int]:
    both_match = 0
    neither_match = 0
    only_rymer_match = 0
    only_kmer_match = 0

    for read in reads:
        for i in range(len(read) - w + 1):
            window = read[i:i+w]
            minimizer = min(window[j:j+k] for j in range(w - k + 1))
            rymer = rymer_transform(minimizer)
            
            rc_minimizer = reverse_complement(minimizer)
            rc_rymer = rymer_transform(rc_minimizer)
            
            minimizer_found = minimizer in minimizer_table or rc_minimizer in minimizer_table
            rymer_found = rymer in rymer_table or rc_rymer in rymer_table

            if minimizer_found and rymer_found:
                both_match += 1
            elif not minimizer_found and not rymer_found:
                neither_match += 1
            elif rymer_found:
                only_rymer_match += 1
            elif minimizer_found:
                only_kmer_match += 1

    return {
        "Both Match": both_match,
        "Neither Match": neither_match,
        "Only Rymer Match": only_rymer_match,
        "Only Minimizer Match": only_kmer_match
    }

# Load the reference sequence and create hash tables
sequence = parse_fasta("rCRS.fa")
k, w = 8, 10
minimizer_table = create_index_table(sequence, k, w)
rymer_sequence = rymer_transform(sequence)
rymer_table = create_index_table(rymer_sequence, k, w)

# Parse the gzipped fastq file and count the matches
reads = parse_gzipped_fastq("bact.fq.gz")
match_counts = count_matches(reads, k, w, minimizer_table, rymer_table)
print(match_counts)

