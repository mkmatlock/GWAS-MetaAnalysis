from pyCSV import *
import geneUtils
import geneVerifier as geneDB

__DEBUG=1

__geneSet = set([])
__traitDict = {}
__studyByGene = {}
__studyByTrait = {}
__pValues = {}
__studyGenes = {}
__gwasCatalogue = pyCSV()

def getGenesForTrait(trait, pvalue = 0.05):
    geneList = set([])
    for studyId in __studyByTrait[trait]:
        if __pValues[studyId] < pvalue and studyId in __studyGenes:
            for geneSym in __studyGenes[studyId]:
                geneList.add(geneSym)
    return geneList


def getStudyInfo(studyId):
    global __gwasCatalogue
    geneString = __gwasCatalogue.get(studyId, 13)
    
    geneList = []
    trait = __gwasCatalogue.get(studyId, 7)
    pvalue = float(__gwasCatalogue.get(studyId, 27))
    
    geneItems = geneString.split(",")
    for geneSym in geneItems:
        geneList.append(geneSym.strip())
    
    return trait, geneList, pvalue

def getGenePValues(geneId):
    global __gwasCatalogue,__studyByGene
    studies = __studyByGene[geneId]
    pvalues = []
    for studyId in studies:
        pvalue = float(__gwasCatalogue.get(studyId, 27))
        pvalues.append(pvalue)
    return pvalues
    
def getPValues(geneId, trait):
    global __gwasCatalogue,__studyByGene,__studyByTrait
    studies = __studyByGene[geneId] & __studyByTrait[trait]
    pvalues = []
    for studyId in studies:
        pvalue = float(__gwasCatalogue.get(studyId, 27))
        pvalues.append(pvalue)
    return pvalues

def __addGene(i, geneSym, geneTrait, __ENABLE_GENE_VERIFICATION, __ENABLE_GENE_UPDATES, updatedGeneSet, invalidGeneSet):
    if __ENABLE_GENE_UPDATES:
        parent = geneDB.findUpdatedSymbol(geneSym)
        if parent != None:
            if __DEBUG > 1:
                print "Replacing",geneSym,"with",parent
            updatedGeneSet.add(geneSym)
            geneSym = parent
    
    if not __ENABLE_GENE_VERIFICATION or geneDB.isApproved(geneSym):
        __geneSet.add(geneSym)
        
        try:
            __studyGenes[i].add(geneSym)
        except KeyError:
            __studyGenes[i] = set([geneSym])
        
        try:
            __traitDict[geneSym].add(geneTrait)
        except KeyError:
            __traitDict[geneSym] = set([geneTrait])
            
        try:
            __studyByGene[geneSym].add(i)
        except KeyError:
            __studyByGene[geneSym] = set([i])
    elif __ENABLE_GENE_VERIFICATION:
       invalidGeneSet.add(geneSym)
    
def init(filename, __ENABLE_GENE_VERIFICATION = 0, __ENABLE_GENE_UPDATES = 0, __INCLUDE_MAPPED_GENES = 0, trait_exclude_file = 0, pfilter = 0.05):
    global __DEBUG, __pValues, __gwasCatalogue, __studyByTrait, __geneSet, __traitDict, __studyByGene
    
    exclude_traits = set([])
    if trait_exclude_file != 0:
        ifile = open(trait_exclude_file,'r')
        for line in ifile:
            exclude_traits.add(line.strip())
        ifile.close()
    
    
    __gwasCatalogue.load(filename, "\t")

    invalidGeneSet = set([])
    updatedGeneSet = set([])
    
    for i in xrange(1, __gwasCatalogue.rows+1):
        geneString = __gwasCatalogue.get(i, 13)
        geneTrait = __gwasCatalogue.get(i, 7)
        
        pvalueText = __gwasCatalogue.get(i, 27)
        pvalue = 0
        
        try:
            pvalue = float(pvalueText)
        except ValueError:
            pvalue = -1
        
        if pvalue > pfilter:
            continue
        if geneTrait in exclude_traits:
            continue
        if geneString==None:
            continue
        if geneString == "":
            continue
        
        __pValues[i] = pvalue
        
        try:
            __studyByTrait[geneTrait].add(i)
        except KeyError:
            __studyByTrait[geneTrait] = set([i])
        
        geneItems = geneString.split(",")
        for item in geneItems:
            
            geneSymbols = item.split(" - ")
            
            for geneSym in geneSymbols:
                geneSym = geneUtils.formatGeneSymbol(geneSym.strip())
                __addGene(i,geneSym, geneTrait, __ENABLE_GENE_VERIFICATION, __ENABLE_GENE_UPDATES, updatedGeneSet, invalidGeneSet)
        
        if __INCLUDE_MAPPED_GENES:
            mappedGenes = __gwasCatalogue.get(i, 14)
            
            mappedItems = mappedGenes.split(";")
            
            for item in mappedItems:
                
                geneSymbols = item.split(" - ")
                
                for geneSym in geneSymbols:
                    geneSym = geneUtils.formatGeneSymbol(geneSym)
                    __addGene(i,geneSym, geneTrait, __ENABLE_GENE_VERIFICATION, __ENABLE_GENE_UPDATES, updatedGeneSet, invalidGeneSet)
        
    invalid_file = open("log\\invalid_gwas.txt",'w')
    
    for geneSym in invalidGeneSet:
        invalid_file.write(geneSym+"\n")
    
    invalid_file.close()
    
    
    if __DEBUG > 0:
        print "\n---------------------------------"
        print "GWAS Invalid Gene Symbols:  ", len(invalidGeneSet)
        print "GWAS Updated Gene Symbols:  ", len(updatedGeneSet)
        print "GWAS Total Genes Remaining: ", len(__geneSet)
        print "---------------------------------\n"