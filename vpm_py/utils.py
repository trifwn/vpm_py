import numpy as np

def divide_processors(
    num_processors: int,
):
    """
    Divide n processors into p groups as evenly as possible
    """
    # Divide processors into NBI, NBJ, NBK so that NBI * NBJ * NBK = number of processors
    # Find the factors of the number of processors
    factors = np.unique(np.array([i for i in range(1, num_processors + 1) if num_processors % i == 0]))
    # Find the subsets of factors that multiply to the number of processors
    subsets = []
    for i in range(len(factors)):
        for j in range(i, len(factors)):
            for k in range(j, len(factors)):
                if factors[i] * factors[j] * factors[k] == num_processors:
                    subsets.append((factors[i], factors[j], factors[k]))
    # Find the subset that has the smallest sum
    min_sum = np.inf
    NBI = 0
    NBJ = 0
    NBK = 0
    for subset in subsets:
        if sum(subset) < min_sum:
            min_sum = sum(subset)
            NBI, NBJ, NBK = subset

    return NBI, NBJ, NBK