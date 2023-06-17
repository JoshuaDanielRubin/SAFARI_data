import random


class ReadGenerator:
    def __init__(self, S, delta):
        self.S = S
        self.delta = delta

    def generate_reads(self, N, L):
        reads = []
        for _ in range(N):
            start = random.randint(0, len(self.S) - 1)
            read = [self.S[i % len(self.S)] for i in range(start, start + L)] # Circular!
            deaminated_bases = set()
            for i in range(len(read)):
                if (read[i] == 'C' or read[i] == 'G') and random.random() < self.delta:
                    if read[i] == 'C':
                        read[i] = 'T'
                    elif read[i] == 'G':
                        read[i] = 'A'
                    deaminated_bases.add(i)
            reads.append((''.join(read), start, deaminated_bases))
        return reads

