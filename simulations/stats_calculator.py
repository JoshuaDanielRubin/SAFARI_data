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
        incorrect_rescue_rate, correct_rescue_rate = self.calculate_rescue_rates()

        return total_minimizer, total_rymer, minimizer_prec, rymer_prec, incorrect_rescue_rate, correct_rescue_rate

    def compute_precision(self, seeds):
        true_positive = sum(1 for seed in seeds if seed[1] is not None)
        false_positive = sum(1 for seed in seeds if seed[1] is None)
        return true_positive / (true_positive + false_positive) if true_positive + false_positive != 0 else 0

    def calculate_rescue_rates(self):
        rymer_incorrect_rescue_count = 0
        rymer_correct_rescue_count = 0
        for seed in self.rymer_seeds:
            pos, correct = seed
            if pos in self.deaminated_bases:  # If the seed maps to a deaminated base
                if not correct:  # If the seed does not map correctly
                    rymer_incorrect_rescue_count += 1
            if correct:  # If the seed maps correctly
                rymer_correct_rescue_count += 1
        rymer_incorrect_rescue_rate = rymer_incorrect_rescue_count / len(self.rymer_seeds) if len(self.rymer_seeds) > 0 else 0
        rymer_correct_rescue_rate = rymer_correct_rescue_count / len(self.rymer_seeds) if len(self.rymer_seeds) > 0 else 0
        return rymer_incorrect_rescue_rate, rymer_correct_rescue_rate

