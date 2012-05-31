
def sort(test_results):
    return sorted(test_results, key = lambda pair: pair[1])

def prune(significance, test_results):
    return [pair for pair in test_results if pair[1] <= significance]

def benjamini(significance, test_results):
    sorted_results = sort(test_results)
    x = len(sorted_results)
    result = [(pair[0], pair[1] * (x - i)) for i, pair in enumerate(sorted_results)]
    return prune(significance, result)

def falseDiscoveryRate(significance, test_results, permutation_results):
    t_categories = set([pair[0] for pair in prune(significance, test_results)])
    f_categories = set([pair[0] for pair in prune(significance, permutation_results)])
    
    return float(len(f_categories)) / float(len(t_categories))

if __name__ == "__main__":
    sig = 0.005
    testset = [('a', 0.003),('b', 0.0006),('c', 0.001),('d', 0.01),('e', 0.0004)] 
    pruned = prune(sig, testset)
    
    print sort(pruned)
    print benjamini(sig, testset)

    print falseDiscoveryRate(sig, testset, [('a', 0.001), ('b', 0.04), ('c', 0.1), ('d', 0.76), ('e', 0.002)])
