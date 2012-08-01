import geneVerifier as geneDB
from pyFASTA import *
from geneUtils import *
from pyCSV import *
import geneUtils
import os

__DEBUG = 1

__targetCatalogue = pyCSV()
__drugCatalogue = pyCSV()
__geneSet = set([])
__geneNames = {}
__drugDict = {}
__drugs = {}

def getDrugsTargetingProteinSet(proteins):
    drugSet = set([])
    for protein in proteins:
        if protein in __drugDict:
            drugSet |= __drugDict[protein]
    return drugSet

def initDruglist(drug_file):
    global __DEBUG
    
    __drugCatalogue.load(drug_file)
    
    for r in xrange(1, __drugCatalogue.rows+1):
        db_id = __drugCatalogue.get(r, 0)
        
        db_name = __drugCatalogue.get(r, 1)
        
        wiki_link = __drugCatalogue.get(r, 17)
        drug_link = __drugCatalogue.get(r, 18)
        
        link = drug_link
        if link == None or link.strip() == "":
            link = wiki_link
        if link.strip() == "":
            link = None
        
        __drugs[db_id] = {'name':db_name, 'link':link}

def initTargets(targets_file, protein_file,__ENABLE_GENE_VERIFICATION=0, __ENABLE_GENE_UPDATES=0):
    global __DEBUG, __targetCatalogue, __geneSet, __geneNames, __drugDict
    
    __targetCatalogue.load(targets_file)
    
    rejectedSet = set([])
    updatedSet = set([])
    
    for r in xrange(1, __targetCatalogue.rows+1):
        geneId = int(__targetCatalogue.get(r,0))
        geneName = geneUtils.formatGeneSymbol(__targetCatalogue.get(r,2))
        
        if geneName != None:
            if __ENABLE_GENE_UPDATES:
                parentSym = geneDB.findUpdatedSymbol(geneName)
                if parentSym != None:
                    updatedSet.add(geneName)
                    geneName = parentSym
                
            if __ENABLE_GENE_VERIFICATION and not geneDB.isApproved(geneName):
                if __DEBUG>2 and geneName != "" or geneName == "papc":
                    print "Rejected:", geneName
                rejectedSet.add(geneName)
                continue
                
            __geneNames[geneId] = geneName
            __geneSet.add(geneName)
            __drugDict[geneName] = set([])
    
    invalid_file = open(os.sep.join(["results","log","invalid_drugbank.txt"]),'w')
    for geneName in rejectedSet:
        invalid_file.write(geneName+"\n")
    invalid_file.close()            
    
    proteins = parseFASTA(protein_file)
    
    __drugSet = set([])
    empty_gene_drug_targets = 0
    for fasta in proteins:
        items = fasta[1].split()
        geneId = int(items[0])
        
        if geneId in __geneNames:
            parenthetical = fasta[1][fasta[1].rfind("(")+1 : fasta[1].rfind(")")]
            
            drugs = parenthetical.split(";")
            
            for drug in drugs:
                drugbankid = drug.strip()
                __drugDict[__geneNames[geneId]].add(drugbankid)
                __drugSet.add(drugbankid)
    
    removable = set([])
    for gene in __geneSet:
        if gene not in __drugDict or len(__drugDict[gene]) == 0:
            removable.add(gene)

            empty_gene_drug_targets += 1
    __geneSet -= removable
    
    removable_drugs = set([])
    for drugbankid in __drugs:
        if drugbankid not in __drugSet:
            removable_drugs.add(drugbankid)
    
    for drugbankid in __drugSet:
        if drugbankid not in __drugs:
            __drugs[drugbankid] = {'name':drugbankid}

    for drugbankid in removable_drugs:
        del __drugs[drugbankid]

    if __DEBUG>0:
        print "\n------------------------------------------"
        print "Invalid Drug Target Gene Symbols:   ", len(rejectedSet)
        print "Updated Drug Target Gene Symbols:   ", len(updatedSet)
        print "Remaining Drug Target Gene Symbols: ", len(__geneSet)
        print "Drugs with targets:                 ", len(__drugSet), len(__drugs)
        print "Removed:", empty_gene_drug_targets, "genes without targeting drugs"
        print "------------------------------------------\n"
       
if __name__ == "__main__":
    print "\nInitializing HGNC Database..."
    geneDB.init(os.sep.join(["data","hgnc","hgnc_symbols.txt"]),1)
    
    print "\nLoading DrugBank drug catalogue..."
    initDruglist(os.sep.join(["data","drugbank","drug_links.csv"]))
    print "\nLoading DrugBank drug target catalogue..."
    initTargets(os.sep.join(["data","drugbank","target_links.csv"]),
    os.sep.join(["data","drugbank","all_target_protein.fasta"]),1,1)

    print "Total Targets:", len(__geneSet)
    print "Total Drugs:", len(__drugs)
    
