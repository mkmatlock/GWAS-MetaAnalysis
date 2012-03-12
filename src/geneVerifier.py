from pyCSV import *
import geneUtils

__DEBUG=0

def __strToList(strItem):
    if strItem == None:
        return []
        
    strItem = strItem.strip()
    if strItem == "":
        return []
    
    items = strItem.split(",")
    
    return [item.strip() for item in items]
    
def __returnIfValid(strItem):
    if strItem == None:
        return None
    return strItem.strip()

__geneDB = 0
__sym_rows = {}

__symbols = set([])
__approved_symbols = set([])
__original_names = {}

# gene DB column names
__hgnc_id = 0
__names = 2
__status = 3
__past_symbols = 4
__synonyms = 5
__chromosome = 6
__accession = 7
__refseq = 8


__gene_symbol_parents = {}
__gene_symbol_synonyms = {}

def findUpdatedSymbol(geneSym):
    global __gene_symbol_parents
    
    geneSym = geneSym.lower()
    if geneSym in __gene_symbol_parents:
        if __DEBUG > 0 and len(__gene_symbol_parents[geneSym])>1:
            print "Gene:",geneSym,"has multiple updated gene names"
            if __DEBUG > 1:
                print __gene_symbol_parents[geneSym]
        
        parentSym = __gene_symbol_parents[geneSym].pop()
        __gene_symbol_parents[geneSym].add(parentSym)
        
        return parentSym
    
    if geneSym in __gene_symbol_synonyms:
        if __DEBUG > 1:
            print __gene_symbol_synonyms[geneSym]
        
        if len(__gene_symbol_synonyms[geneSym])>0:
            synSym = __gene_symbol_synonyms[geneSym].pop()
            __gene_symbol_synonyms[geneSym].add(synSym)
        
            return synSym
    
    return None

    
def getHGNCId(geneSym):
    global __geneDB, __sym_rows, __hgnc_id
    return __returnIfValid(__geneDB.get(__sym_rows[geneSym.lower()], __hgnc_id))

    
def getName(geneSym):
    global __geneDB, __sym_rows, __names
    return __returnIfValid(__geneDB.get(__sym_rows[geneSym.lower()], __names))
    
    
def isValid(geneSym):
    global __symbols
    return geneSym.lower() in __symbols
    
    
def isApproved(geneSym):
    global __approved_symbols
    return geneSym.lower() in __approved_symbols    
    
    
def __isApproved(geneSym):
    val = __geneDB.get(__sym_rows[geneSym.lower()], __status)
    return val != None and "approved" == val.strip().lower()
    
    
def getStatus(geneSym):
    global __geneDB, __sym_rows, __status
    return __returnIfValid(__geneDB.get(__sym_rows[geneSym.lower()], __status))

    
def getPastSymbols(geneSym):
    global __geneDB, __sym_rows, __past_symbols
    return [geneUtils.formatGeneSymbol(item) for item in __strToList(__geneDB.get(__sym_rows[geneSym.lower()], __past_symbols))]

    
def getSynonyms(geneSym):
    global __geneDB, __sym_rows, __synonyms
    return [geneUtils.formatGeneSymbol(item) for item in __strToList(__geneDB.get(__sym_rows[geneSym.lower()], __synonyms))]
    
    
def getChromosome(geneSym):
    global __geneDB, __sym_rows, __accession
    return __returnIfValid(__geneDB.get(__sym_rows[geneSym.lower()], __accession))
    
    
def getAccessionNumbers(geneSym):
    global __geneDB, __sym_rows, __accession
    return __strToList(__geneDB.get(__sym_rows[geneSym.lower()], __accession))  

    
def getRefSeqIds(geneSym):
    global __geneDB, __sym_rows, __refseq
    return __strToList(__geneDB.get(__sym_rows[geneSym.lower()], __refseq))

    
def init(filename, __EXCLUDE_PSEUDOGENES=0):
    global __geneDB, __symbols, __gene_symbol_parents, __gene_symbol_synonyms
    
    __geneDB = pyCSV()
    __geneDB.load(filename,"\t")
    
    for row in xrange(1, __geneDB.rows+1):
        
        if __EXCLUDE_PSEUDOGENES and __geneDB.get(row, 2).count("pseudogene")>0:
            continue
        
        symbol = geneUtils.formatGeneSymbol(__geneDB.get(row, 1))
        
        __original_names[symbol] = __geneDB.get(row, 1)
        __symbols.add(symbol)
        
        __sym_rows[symbol] = row
        
        if __isApproved(symbol):
            __approved_symbols.add(symbol)
        
    print "Loaded:", len(__approved_symbols), "approved gene symbols..."
        
    for symbol in __approved_symbols:
        past_symbols = getPastSymbols(symbol)
        
        for child in past_symbols:
            try:
                __gene_symbol_parents[child].add(symbol)
            except KeyError:
                __gene_symbol_parents[child] = set([symbol])
    
        for synSym in getSynonyms(symbol):
            try:
                __gene_symbol_synonyms[synSym].add(symbol)
            except KeyError:
                __gene_symbol_synonyms[synSym] = set([symbol])
                
    for symbol in __gene_symbol_synonyms:
        remove = []
        
        for synGene in __gene_symbol_synonyms[symbol]:
            if synGene not in __approved_symbols:
                remove.append(synGene)
        
        for r in remove:
            __gene_symbol_synonyms[symbol].remove(r)