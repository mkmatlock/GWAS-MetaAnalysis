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
    return fisher.compute(a,b,c,d)

def simultaneousPermutation(n, category_sizes, testset, iterations):
    avg_p_values = [0] * len(category_sizes)
    
    print "Running permutation test on %d categories drawn from %d items for %d iterations" % (len(category_sizes), n, iterations)
    pbar = ProgressBar()
    pbar.setMaximum(iterations)
    pbar.updateProgress(0)

    sampleset = xrange(0, n)
    for i in xrange(0, iterations):
        if i % 5 == 0:
            pbar.updateProgress(i)
        for j, c in enumerate(category_sizes):
            s1 = set(random.sample(sampleset, c))
            avg_p_values[j] += computeOverlapPval(s1, testset, n)

    pbar.finalize()

    return [v / float(iterations) for v in avg_p_values]

def simultaneousPermutationIncludeTestSet(n, category_sizes, testset_size, iterations):
    avg_p_values = [0] * len(category_sizes)

    sampleset = xrange(0, n)
    for i in xrange(0, iterations):
        for j, c in enumerate(category_sizes):
            s1 = set(random.sample(sampleset, c))
            s2 = set(random.sample(sampleset, testset_size))
            avg_p_values[j] += computeOverlapPval(s1, s2, n)

    return [v / float(iterations) for v in avg_p_values]

if __name__ == "__main__":
    n = 50
    iterations = 1000

    fisher.init(n)

    testset = set([0, 1, 4, 7, 2, 10, 40, 32])

    categories = [set([0,1]),set([1,4,7,0,2,10,40]),set([10,2,40]),set([32])]
    category_sizes = [len(s) for s in categories]

    cat_pvalues = [computeOverlapPval(s, testset, n) for s in categories]
    perm_pvalues = simultaneousPermutation(n, category_sizes, testset, iterations)
    perm_pvalues2 = simultaneousPermutationIncludeTestSet(n, category_sizes, len(testset), iterations)

    for i in xrange(0, len(categories)):
        print "%s || %.7f     %.7f     %.7f" % ("%d" % (i), cat_pvalues[i], perm_pvalues[i], perm_pvalues2[i])
