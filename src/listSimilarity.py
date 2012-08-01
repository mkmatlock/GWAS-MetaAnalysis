import os,sys


listA = []
listB = []


for i in xrange(1, len(sys.argv)):
    if sys.argv[i] == "--appendA":
        fileA = open(sys.argv[i+1],'r')
        for line in fileA:
            item = line.strip()
            if item != "":
                listA.append(item)
        fileA.close()
        
    if sys.argv[i] == "--appendB":
        fileB = open(sys.argv[i+1],'r')
        for line in fileB:
            item = line.strip()
            if item != "":
                listB.append(item)
        fileB.close()

if len(listA) == 0 or len(listB) == 0:
    print "Nothing to do..."
    sys.exit(1)
        
def isdigit(a):
    return a == '0' or a == '1' or a == '2' or a == '3' or a == '4' or a == '5' or a =='6' or a == '7' or a == '8' or a == '9'


def damerau_levenshtein_distance(s1, s2):
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)
    for i in xrange(-1,lenstr1+1):
        d[(i,-1)] = i+1
    for j in xrange(-1,lenstr2+1):
        d[(-1,j)] = j+1
 
    for i in xrange(lenstr1):
        for j in xrange(lenstr2):
            if s1[i] == s2[j]:
                cost = 0
            elif isdigit(s1[i]) and s1[i] != s2[j]:
                cost = 3
            elif i < lenstr1-1 and j < lenstr2-1 and s1[i] == s2[j+1] and s2[j] == s1[i+1]:
                cost = 0
            else:
                cost = 1
            d[(i,j)] = min(
                           d[(i-1,j)] + 2, # deletion
                           d[(i,j-1)] + 2, # insertion
                           d[(i-1,j-1)] + cost, # substitution
                          )
            if i and j and s1[i]==s2[j-1] and s1[i-1] == s2[j]:
                d[(i,j)] = min (d[(i,j)], d[i-2,j-2] + cost) # transposition
 
    return d[lenstr1-1,lenstr2-1]

def similarity(A,B):
    if A == B:
        return 0
    
    if B.startswith(A) or A.startswith(B):
        return 0.5
    
    return damerau_levenshtein_distance(A,B)

pairs = []
for a in listA:
    for b in listB:
        pairs.append((a,b,similarity(a,b)))
        
sorted_pairs = sorted(pairs, key=lambda item: item[2])

lval = 0
for item in sorted_pairs:
    if item[2] != lval:
        print "---------------------------------------------"
        lval = item[2]
    print "%-20s %-20s %.1f" % tuple(item)