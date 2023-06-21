class StatsCalculator:
    def __init__(self, minimizer_seeds, rymer_seeds, deaminated_bases, k):
        self.minimizer_seeds = minimizer_seeds
        self.rymer_seeds = rymer_seeds
        self.deaminated_bases = deaminated_bases
        self.k = k

    def compute_stats(self):
        total_minimizer = len(self.minimizer_seeds)
        total_rymer = len(self.rymer_seeds)
        
        minimizer_prec = self.compute_precision(self.minimizer_seeds)
        rymer_prec = self.compute_precision(self.rymer_seeds)
        rymer_spuriousness = self.compute_rymer_spuriousness()
        rymer_recovery = self.compute_rymer_recovery()

        return total_minimizer, total_rymer, minimizer_prec, rymer_prec, rymer_spuriousness, rymer_recovery

    def compute_precision(self, seeds):
        true_positive = sum(1 for seed in seeds if seed[1] is not None)
        false_positive = sum(1 for seed in seeds if seed[1] is None)
        return true_positive / (true_positive + false_positive) if true_positive + false_positive != 0 else 0

    def compute_rymer_spuriousness(self):
        false_positive = sum(1 for seed in self.rymer_seeds if seed[1] is None)
        return false_positive / len(self.rymer_seeds) if len(self.rymer_seeds) != 0 else 0

    def compute_rymer_recovery(self):
        true_positive = sum(1 for seed in self.rymer_seeds if seed[1] is not None)
        return true_positive / len(self.deaminated_bases) if len(self.deaminated_bases) != 0 else 0

