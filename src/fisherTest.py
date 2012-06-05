import numpy as np
import math
import sys

__primeFactors = {}
__primeSeive = []
__primeSet = set([])

def factorize(n):
    primes = {}
    for i in xrange(0, len(__primeSeive)):
        p = __primeSeive[i]

        if p > n:
            break
        
        while n % p == 0:
            try:
                primes[p] += 1
            except KeyError:
                primes[p] = 1
            n = n / p
    
    return primes

def copyFactors(n, factors):
    global __primeFactors
    if n == 0 or n == 1:
        return
    for p in __primeFactors[n]:
        try:
            factors[p] += __primeFactors[n][p]
        except KeyError:
            factors[p] = __primeFactors[n][p]

def computeFactorsUpTo(n):
    global __primeFactors, __primeSet
    __primeFactors = {}

    for i in xrange(2,n+1):
        factors = {}
        if i in __primeSet:
            factors[i] = 1
        else:
            nRt = int(math.sqrt(i))
            for j in xrange(2, nRt+1):
                if i % j == 0:
                    copyFactors(j, factors)
                    copyFactors(i/j, factors)
                    __primeFactors[i] = factors
                    break
        __primeFactors[i] = factors


def computeSeive(n):
    global __primeSeive
    queue = range(2,n)
    __primeSeive = []

    while len(queue) > 0:
        p = queue.pop(0)
        __primeSeive.append(p)
        __primeSet.add(p)
        
        i = 2
        while (i*p) < n:
            try:
                queue.remove(i*p)
            except ValueError:
                pass
            i+=1

def printPrimeFactorization(n):
    factors = __primeFactors[n]

    keysort = sorted(factors.keys())

    expression = str(n) + " = "

    for p in keysort:
        n = factors[p]
        if n == 1:
            expression += "("+str(p)+")"
        else:
            expression += "("+str(p)+"^"+str(n)+")"
    print expression

def eliminateFactors(factors1, factors2):
    keySet1 = set(factors1.keys())
    keySet2 = set(factors2.keys())

    for p in keySet1 & keySet2:
        n1 = factors1[p]
        n2 = factors2[p]

        if n1 > n2:
            del factors2[p]
            factors1[p] = n1 - n2
        elif n1 < n2:
            del factors1[p]
            factors2[p] = n2 - n1
        else:
            del factors1[p]
            del factors2[p]

def computeV1(a,b,c,d):
    def addAllFactorsBetween((n,m),factors):
        for i in range(n,m+1):
            copyFactors(i, factors)

    n = a+b+c+d
    denom_factorials = [(a+c+1,n),(2,d)]
    numer_factorials = [(a+1,a+b),(c+1,c+d),(b+1,b+d)]

    numer_factors = {}
    denom_factors = {}
    
    [addAllFactorsBetween(R, denom_factors) for R in denom_factorials]
    [addAllFactorsBetween(R, numer_factors) for R in numer_factorials]
    
    eliminateFactors(numer_factors, denom_factors)
    
    numer = sorted([math.pow(p,numer_factors[p]) for p in numer_factors])
    denom = sorted([math.pow(p,denom_factors[p]) for p in denom_factors])
    numer.reverse()
    denom.reverse()
    
    return __compute_quotient(numer, denom)

def computeV2(a,b,c,d):
    def addAllItemsBetween((n,m),factors):
        for i in range(n,m+1):
            try:
                factors[i] += 1
            except KeyError:
                factors[i] = 1

    n = a+b+c+d
    denom_factorials = [(a+c+1,n),(2,d)]
    numer_factorials = [(a+1,a+b),(c+1,c+d),(b+1,b+d)]

    numer_factors = {}
    denom_factors = {}
    
    [addAllItemsBetween(R, denom_factors) for R in denom_factorials]
    [addAllItemsBetween(R, numer_factors) for R in numer_factorials]
    
    eliminateFactors(numer_factors, denom_factors)
    
    numer = [math.pow(p,numer_factors[p]) for p in numer_factors]
    denom = [math.pow(p,denom_factors[p]) for p in denom_factors]
    
    return __compute_quotient(numer, denom)


def __compute_quotient(numer, denom):
    p=1.0
    while len(numer) > 0 and len(denom) > 0:
        if p>1.0:
            p /= float(denom[0])
            denom.pop(0)
        else:
            p *= float(numer[0])
            numer.pop(0)
        
    while len(numer) > 0:
        p *= float(numer[0])
        numer.pop(0)
    
    while len(denom) > 0:
        p /= float(denom[0])
        denom.pop(0)
    
    return p

def reduceTable(a,b,c,d):
    while a > 0:
        a = a-1
        b = b+1
        c = c+1
        d = d-1
    return a,b,c,d

def restoreTable(i, a, b, c, d, fisher):
    while i > 0:
        fisher *= float(b * c) / float((a+1)*(d+1))
        a+=1
        b-=1
        c-=1
        d+=1
        
        i-=1

    return fisher

__cached = {}

def compute(a,b,c,d):
    ar,br,cr,dr = reduceTable(a,b,c,d)
    
    key = (ar,br,cr,dr)
    if key in __cached:
        fisher = __cached[key]
    else:
        fisher = computeV2(ar,br,cr,dr)
        __cached[key] = fisher

    return restoreTable(a, ar,br,cr,dr, fisher)


def init(n):
    computeSeive(n)
    computeFactorsUpTo(n)

if __name__ == "__main__":
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    c = int(sys.argv[3])
    d = int(sys.argv[4])
    n = a+b+c+d
    
    print "computing factors to %d..." % (n)
    computeSeive(n)
    computeFactorsUpTo(n)
    
    print "compute fisher exact test a=%d, b=%d, c=%d, d=%d" % (a,b,c,d)
    print compute(a,b,c,d)

    import time

    print "scale test"
    start = time.clock()
    for i in xrange(0,100):
        compute(a+i,b-i,c-i,d+i)
    end = time.clock()
    print "time elapsed:", (end-start)
