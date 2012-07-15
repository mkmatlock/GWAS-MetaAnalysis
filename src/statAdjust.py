
def sort(test_results):
    return sorted(test_results, key = lambda pair: pair[1])

def prune(significance, test_results):
    return [pair for pair in test_results if pair[1] <= significance]

def bonferoni(significance, test_results):
    sorted_results = sort(test_results)
    x = len(sorted_results)
    return [(pair[0], pair[1] * (x - i)) for i, pair in enumerate(sorted_results)]


def falseDiscoveryRate(significance_levels, permutation_results, N_h):
    return [ __internal_falseDiscoveryRate(significance, permutation_results, N_h) for significance in significance_levels ]

def __internal_falseDiscoveryRate(significance, permutation_results, N_h):
    f_categories = set([pair[0] for pair in prune(significance, permutation_results)])

    return float(len(f_categories)) / float(N_h)

if __name__ == "__main__":
    sig = 0.005
    testset = [('a', 0.003),('b', 0.0006),('c', 0.001),('d', 0.01),('e', 0.0004)] 
    pruned = prune(sig, testset)

    print sort(pruned)
    print bonferoni(sig, testset)

    print falseDiscoveryRate([0.05, sig, 0.0005], [('a', 0.001), ('b', 0.04), ('c', 0.1), ('d', 0.76), ('e', 0.002)], len(testset))
