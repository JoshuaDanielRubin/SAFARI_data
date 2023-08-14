from collections import defaultdict

class UniqueMinimizerSketch:
    def __init__(self, S, W, k):
        self.S = S
        self.W = W
        self.k = k
        self.minimizer_index, self.rymer_index, self.rymer_minimizer_map = self.create_minimizer_index()

    def create_minimizer_index(self):
        minimizer_index = defaultdict(list)
        rymer_index = defaultdict(list)
        rymer_minimizer_map = defaultdict(set)

        # First pass: construct the full minimizer index.
        for i in range(len(self.S) - self.W + 1):
            window = self.S[i:i + self.W]
            for j in range(len(window) - self.k + 1):
                kmer = window[j:j + self.k]
                minimizer = min(kmer, self.reverse_complement(kmer))
                minimizer_index[minimizer].append(i + j)

        # Second pass: select the most unique minimizer for each window.
        unique_minimizer_index = defaultdict(list)
        for i in range(len(self.S) - self.W + 1):
            window = self.S[i:i + self.W]
            minimizers = []
            for j in range(len(window) - self.k + 1):
                kmer = window[j:j + self.k]
                minimizer = min(kmer, self.reverse_complement(kmer))
                minimizers.append((minimizer, i + j))

            # Select the most unique minimizer in the window.
            unique_minimizer = min(minimizers, key=lambda x: (len(minimizer_index[x[0]]), x[0]))
            minimizer, pos = unique_minimizer
            unique_minimizer_index[minimizer].append(pos)

            # Create rymer and update rymer index and rymer minimizer map.
            rymer = ''.join(['R' if base in 'GA' else 'Y' for base in minimizer])
            rymer_index[rymer].append(pos)
            rymer_minimizer_map[rymer].add(minimizer)

        return unique_minimizer_index, rymer_index, rymer_minimizer_map

    @staticmethod
    def reverse_complement(kmer):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(complement[base] for base in reversed(kmer))

    def find_seeds(self, index, read, origin):
        seeds = []
        for i in range(len(read) - self.k + 1):
            kmer = read[i:i + self.k]
            if kmer in index:
                seeds.extend([(pos, pos == origin + i) for pos in index[kmer]])
        return seeds


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
            minimizers = []
            for j in range(len(window) - self.k + 1):
                kmer = window[j:j + self.k]
                minimizer = min(kmer, self.reverse_complement(kmer))
                minimizers.append((minimizer, i + j))

            # select the minimizer in the window
            minimizer, pos = min(minimizers)

            # update the minimizer index
            minimizer_index[minimizer].append(pos)

            # create rymer and update rymer index and rymer minimizer map
            rymer = ''.join(['R' if base in 'GA' else 'Y' for base in minimizer])
            rymer_index[rymer].append(pos)
            if rymer not in rymer_minimizer_map:
                rymer_minimizer_map[rymer] = set()
            rymer_minimizer_map[rymer].add(minimizer)

        return minimizer_index, rymer_index, rymer_minimizer_map

    @staticmethod
    def reverse_complement(kmer):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        return "".join(complement[base] for base in reversed(kmer))

    def find_seeds(self, index, read, origin):
        seeds = []
        for i in range(len(read) - self.k + 1):
            kmer = read[i:i + self.k]
            if kmer in index:
                seeds.extend([(pos, pos == origin + i) for pos in index[kmer]])
        return seeds

