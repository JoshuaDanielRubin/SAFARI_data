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
        correct_rescue_rate = self.calculate_correct_rescue_rate()

        return total_minimizer, total_rymer, minimizer_prec, rymer_prec, correct_rescue_rate

    def compute_precision(self, seeds):
        true_positive = sum(1 for seed in seeds if seed[1] is not None)
        false_positive = sum(1 for seed in seeds if seed[1] is None)
        return true_positive / (true_positive + false_positive) if true_positive + false_positive != 0 else 0

    def calculate_correct_rescue_rate(self):
        rymer_correct_rescue_count = sum(1 for seed in self.rymer_seeds if seed[1])
        return rymer_correct_rescue_count / len(self.rymer_seeds) if len(self.rymer_seeds) > 0 else 0
