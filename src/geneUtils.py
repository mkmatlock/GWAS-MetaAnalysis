import geneVerifier as geneDB
import math

__DEBUG=0


def contingentChiSquare(a,b,c,d):
    return math.pow(a*d-b*c,2) * (a + b + c + d) / ((a+b)*(c+d)*(b+d)*(a+c))
    
replacement_dict = {"hat":"hat1", "cdkal":"cdkal1"}

def geneFrequency(geneFreqDict, geneSet):
    for gene in geneSet:
        try:
            geneFreqDict[gene]+=1
        except KeyError:
            geneFreqDict[gene] = 1

def addFilterFrequency(geneFreq, __cutoff=1):
    geneSet = set([])
    duplicates = 0
    for gene in geneFreq:
        if geneFreq[gene] >= __cutoff:
            geneSet.add(gene)
            duplicates += geneFreq[gene] - 1
    return duplicates, geneSet
    

def formatGeneSymbol(geneSym):
    geneSym=geneSym.replace("-","").replace("/","_")
    
    tilde = geneSym.find("~")
    star  = geneSym.find("*")
    
    if tilde > -1:
        geneSym = geneSym[0:tilde]
    
    if star > -1:
        geneSym = geneSym[0:star]
    
    geneSym=geneSym.lower()
    
    if geneSym in replacement_dict:
        return replacement_dict[geneSym]
    
    return geneSym

def mergeColumns(csv, c1, c2, r1=0, r2 = -1):
    global __DEBUG
    if r2 == -1:
        r2 = csv.rows+1
    list = []
    conflicts=[]
    for r in xrange(r1,r2):
        v1 = csv.get(r,c1)
        v2 = csv.get(r,c2)
        
        if v1==v2 and v1=="" or v1==None:
            pass
        elif v1==v2:
            list.append(v1)
        elif v1=="" or v1==None:
            list.append(v2)
        elif v2=="" or v2==None:
            list.append(v1)
        elif v1!=v2:
            if __DEBUG>2:
                print "WARNING: columns [%d,%d] do not match at row %d with values: \"%s\" != \"%s\"" % (c1,c2,r,v1,v2)
            conflicts.append((v1,v2))
            
            
    
    return list, conflicts
    
def columnToList(csvfile, col, start = 0, end = -1):
    if end == -1:
        end = csvfile.rows+1
    list = []
    
    for row in xrange(start, end):
        val = csvfile.get(row, col)
        
        if val != None and val != "":
            list.append(val)
            
    return list
    
def addAll(set, list):
    dupes=0
    
    for item in list:
        if item in set:
            dupes+=1
            if __DEBUG>3:
                print "WARNING: Duplicate item:", item
        else:
            set.add(item)
    return dupes

def removeInvalidGenes(geneSet):
    remove = []
    for gene in geneSet:
        if not geneDB.isApproved(gene):
            if __DEBUG>2:
                print "Gene:", gene,"is not valid..."
            remove.append(gene)
            
    invalid_file = open("log\\invalid_genelist.txt",'w')
    for r in remove:
        geneSet.remove(r)
        invalid_file.write(r + "\n")
    invalid_file.close()
        
    if __DEBUG>0:
        print "\n----------------------------"
        print "Gene Set Correction Results:"
        print "Removed:     ", len(remove)
        print "----------------------------\n"
    
def updateGeneSet(geneSet):
    global __DEBUG
    # return complete lists
    remove = []
    add = []
    for gene in geneSet:
        # check in the deprecated listings
        parent = geneDB.findUpdatedSymbol(gene)
        if parent != None:
            if __DEBUG>1:
                print "Replacing",gene,"with",parent
            add.append(parent)
            remove.append(gene)
            
    for r in remove:
        geneSet.remove(r)
    for a in add:
        geneSet.add(a)
        
    if __DEBUG>0:
        print "\n----------------------------"
        print "Gene Set Update Results:"
        print "Updated:     ", len(remove)
        print "----------------------------\n"