import sys 
import fisherTest as fisher
import random
from drugbankDatabase import ProgressBar

def __internal_compute(a, c1, c2, n, iterations):
    items = xrange(0,n)
    n_g = 0
    for i in xrange(0, iterations):
        s1 = set(random.sample(items, c1))
        s2 = set(random.sample(items, c2))
        l = len(s1 & s2)
        if a <= l:
            n_g += 1
    return n_g / iterations

def contingencyMatrix(a,b,c,d, iterations):
    n = a+b+c+d
    return __internal_compute(a, a+b, a+c, n, iterations)

def computeOverlapPval(s1, s2, n):
    a = len(s1 & s2)
    b = len(s1 - s2)
    c = len(s2 - s1)
    d = n - (a + b + c)
    probability_value = fisher.compute(a,b,c,d)
    return fisher.significance(probability_value, a,b,c,d)


def simultaneousPermutation(significance_level, n, categories, testset_size, iterations):
    
    print "Running permutation test on %d categories with test set size of %d drawn from %d items for %d iterations" % (len(categories), testset_size, n, iterations)
    pbar = ProgressBar()
    pbar.setMaximum(iterations)
    pbar.updateProgress(0)

    total_sig_hypotheses = 0

    sampleset = xrange(0, n)
    for i in xrange(0, iterations):
        if i % 1 == 0:
            pbar.updateProgress(i)
        test_set = set(random.sample(sampleset, testset_size))
        total_sig_hypotheses += len([j for j, s_c in enumerate(categories) if computeOverlapPval(s_c, test_set, n) < significance_level])

    pbar.finalize()

    return float(total_sig_hypotheses) / float(iterations)

if __name__ == "__main__":
    n = 40
    iterations = 1000

    fisher.init(n)

    testset = set([0, 1, 4, 7, 2, 10, 40, 32])

    categories = [set([0,1]),set([1,4,7,0,2,10,40]),set([10,2,40]),set([32])]

    cat_pvalues = [computeOverlapPval(s, testset, n) for s in categories]
    perm_DR = simultaneousPermutation(0.05, n, categories, len(testset), iterations)

    print "E[P_rand | H0] = %.7f" % (perm_DR)
    for i in xrange(0, len(categories)):
        print "%s || %.7f" % ("%d" % (i), cat_pvalues[i])
