from collections import defaultdict


class MinimizerSketch:
    def __init__(self, S, W, k):
        self.S = S
        self.W = W
        self.k = k
        self.minimizer_index, self.rymer_index, self.rymer_minimizer_map = self.create_minimizer_index()

    def create_minimizer_index(self):
        minimizer_index = defaultdict(list)
        rymer_index = defaultdict(list)
        rymer_minimizer_map = {}

        for i in range(len(self.S) - self.W + 1):
            window = self.S[i:i + self.W]
            for j in range(len(window) - self.k + 1):
                kmer = window[j:j + self.k]
                minimizer = min(kmer, self.reverse_complement(kmer))
                rymer = ''.join(['R' if base in 'GA' else 'Y' for base in kmer])
                minimizer_index[minimizer].append(i + j)
                rymer_index[rymer].append(i + j)
                if rymer not in rymer_minimizer_map:
                    rymer_minimizer_map[rymer] = set()
                rymer_minimizer_map[rymer].add(minimizer)

        return minimizer_index, rymer_index, rymer_minimizer_map

    @staticmethod
    def reverse_complement(kmer):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement[base] for base in reversed(kmer))

    def find_seeds(self, index, read, origin):
        seeds = []
        for i in range(len(read) - self.k + 1):
            kmer = read[i:i + self.k]
            if kmer in index:
                seeds.extend([(pos, pos == origin + i) for pos in index[kmer]])
        return seeds

