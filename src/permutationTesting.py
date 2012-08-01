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

def simultaneousPermutationWithMultipleSignificanceLevels(significance_levels, n, categories, testset_size, iterations):
    print "Running permutation test on %d categories with test set size of %d drawn from %d items for %d iterations" % (len(categories), testset_size, n, iterations)
    pbar = ProgressBar()
    pbar.setMaximum(iterations)
    pbar.updateProgress(0)

    total_sig_hypotheses = [0] * len(significance_levels)
    num_sig_levels = len(significance_levels)

    sampleset = xrange(0, n)
    for i in xrange(0, iterations):
        if i % 1 == 0:
            pbar.updateProgress(i)
        test_set = set(random.sample(sampleset, testset_size))
        for j, s_c in enumerate(categories):
            pval = computeOverlapPval(s_c, test_set, n)

            for k in xrange(0, num_sig_levels):
                if pval <= significance_levels[k]:
                    total_sig_hypotheses[k] += 1

    pbar.finalize()

    return [ float(total_sig) / float(iterations) for total_sig in total_sig_hypotheses ]

def simultaneousPermutation(significance_level, n, categories, testset_size, iterations):
    result = simultaneousPermutationWithMultipleSignificanceLevels([significance_level], n, categories, testset_size, iterations)
    return result[0]


if __name__ == "__main__":
    n = 40
    iterations = 1000

    fisher.init(n)

    testset = set([0, 1, 4, 7, 2, 10, 40, 32])

    categories = [set([0,1]),set([1,4,7,0,2,10,40]),set([10,2,40]),set([32])]

    cat_pvalues = [computeOverlapPval(s, testset, n) for s in categories]
    levels = [0.1, 0.05, 0.01, 0.005]
    perm_DR = simultaneousPermutationWithMultipleSignificanceLevels(levels, n, categories, len(testset), iterations)

    print "FDR:"
    for i in xrange(0,4):
        print "E[V/R] = %.7f at %.4f" % (perm_DR[i], levels[i])

    print "Sample permutation results:"
    for i in xrange(0, len(categories)):
        print "%s || %.7f" % ("%d" % (i), cat_pvalues[i])
