class StatsCalculator:
    def __init__(self, minimizer_seeds, rymer_seeds, deaminated_bases, k):
        self.minimizer_seeds = minimizer_seeds
        self.rymer_seeds = rymer_seeds
        self.deaminated_bases = deaminated_bases
        self.k = k

    def compute_stats(self):
        total_minimizer_seeds, minimizer_precision = self.compute_seed_stats(self.minimizer_seeds)
        total_rymer_seeds, rymer_precision = self.compute_seed_stats(self.rymer_seeds)
        rymer_spuriousness = self.compute_rymer_spuriousness()
        rymer_recovery_rate = self.compute_rymer_recovery()
        return total_minimizer_seeds, total_rymer_seeds, minimizer_precision, rymer_precision, rymer_spuriousness, rymer_recovery_rate

    @staticmethod
    def compute_seed_stats(seeds):
        total_seeds_found = len(seeds)
        correct_seeds = sum(correct for _, correct in seeds)
        if total_seeds_found > 0:
            precision = correct_seeds / total_seeds_found
        else:
            precision = 0
        return total_seeds_found, precision

    def compute_rymer_spuriousness(self):
        minimizer_kmers = {pos for pos, _ in self.minimizer_seeds}
        rymer_kmers = {pos for pos, correct in self.rymer_seeds if correct}

        # Kmers where the minimizer does not match but the rymer does, and is correct
        spurious_rymer_kmers = rymer_kmers.difference(minimizer_kmers)

        # Fraction of these kmers over the total number of rymer seeds
        if rymer_kmers:
            rymer_spuriousness = 1 - (len(spurious_rymer_kmers) / len(rymer_kmers))
        else:
            rymer_spuriousness = 0

        return rymer_spuriousness

    def compute_rymer_recovery(self):
        minimizer_kmers = {pos for pos, _ in self.minimizer_seeds}
        rymer_kmers = {pos for pos, correct in self.rymer_seeds if correct}

        # Kmers where the minimizer does not match but the rymer does, and is correct
        spurious_rymer_kmers = rymer_kmers.difference(minimizer_kmers)

        # Of these, the ones that can be explained by deamination
        deaminated_rymer_kmers = {pos for pos in spurious_rymer_kmers for i in range(self.k) if pos + i in self.deaminated_bases}

        # Fraction of these kmers over the total number of spurious rymer kmers
        if spurious_rymer_kmers:
            rymer_recovery_rate = len(deaminated_rymer_kmers) / len(spurious_rymer_kmers)
        else:
            rymer_recovery_rate = 0

        return rymer_recovery_rate

