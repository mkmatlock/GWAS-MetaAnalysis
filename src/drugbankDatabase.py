import os
import re
import sys
import geneVerifier as geneDB
import geneUtils
from pyCSV import pyCSV
import drugbankCatalog as drugCat

__drugDict = {}
__drugs = {}
__targets = {}
__excluded_targets = 0
__target_names = {}
__targets_unnamed = 0
__geneSet = set([])
START = 0
END = 1
CLOSURE = 2


# private variables
__l_drug = 0
__l_target = 0

__capture = 0
__cap_string = 0


# function resolves
__strstrip = str.strip
__strfind = str.find
__strstartswith = str.startswith
__strendswith = str.endswith
__strsplit = str.split

__listappend = list.append
__listpop = list.pop
__DEBUG=0

def getDrugsTargetingProteinSet(proteins):
    drugs = set([])

    for gene in proteins:
        if gene in __drugDict:
            drugs |= __drugDict[gene]

    return drugs

def mapTargetNames(targets_file,__ENABLE_GENE_UPDATES=1,__ENABLE_GENE_VERIFICATION=1):
    global __targets_unnamed, __targets, __target_names, __geneSet
    
    __targetCatalogue = pyCSV()
    __targetCatalogue.load(targets_file)
    rejectednames = 0

    representedids = set([])

    for r in xrange(1, __targetCatalogue.rows+1):
        targetId   = int(__targetCatalogue.get(r,0))
        targetGene = geneUtils.formatGeneSymbol(__targetCatalogue.get(r,2))
        
        representedids.add(targetId)

        if targetGene != None:
            parentSym = geneDB.findUpdatedSymbol(targetGene)
            if __ENABLE_GENE_UPDATES and parentSym != None:
                if __DEBUG>1:
                    print "Updated:", targetGene, "to:", parentSym
                targetGene = parentSym
            if __ENABLE_GENE_VERIFICATION and not geneDB.isApproved(targetGene):
                if __DEBUG>1:
                    print "Rejected:", targetGene
                rejectednames+=1
                continue

            __geneSet.add(targetGene)
            __target_names[targetId] = targetGene

    
    for targetId in __targets:
        if targetId not in __target_names.keys():
            __targets_unnamed += 1

    if __DEBUG>0:
        print "Rejected:        ", rejectednames
        print "Unrepresented:   ", len(set(__targets.keys()) - representedids)
    
    __drugSet = set([])
    for drugbankid in __drugs:
        drug = __drugs[drugbankid]

        for target in drug['targets']:
            targetId = target['partner']
            if targetId in __target_names:
                targetGene = __target_names[targetId]
                
                __drugSet.add(drugbankid)
                try:
                    __drugDict[targetGene].add(drugbankid)
                except KeyError:
                    __drugDict[targetGene] = set([drugbankid])
    
    removable_drugs = set([])
    for drugbankid in __drugs:
        if drugbankid not in __drugSet:
            removable_drugs.add(drugbankid)

    for drugbankid in removable_drugs:
        del __drugs[drugbankid]

    lenbefore = len(__geneSet)
    __geneSet = __geneSet & set(__drugDict.keys())
    lenafter = len(__geneSet)
    if __DEBUG>0:
        print "Removed", (lenafter - lenbefore), "untargeted gene names"
    
    print "Total drugs with targets:   ", len(__drugSet), len(__drugs)
    print "Total geneset size:         ", len(__geneSet)

def startCapture():
    global __capture, __cap_string
    __capture = 1
    __cap_string = ""

def capture(text):
    global __capture, __cap_string
    if __capture:
        __cap_string += text

def endCapture():
    global __capture, __cap_string
    __capture = 0
    tmp = __cap_string
    __cap_string = 0
    return tmp

def checkParent(tag, stack):
    if len(stack) > 1:
        return tag == stack[-2]
    return False

def handleStart(tag, args, stack):
    global __l_drug, __l_target, __targets
    if tag == 'drug' and checkParent('drugs', stack):
        __l_drug = {}

    elif checkParent('drug', stack):
        if tag == 'drugbank-id' or tag == 'name':
            startCapture()
        if tag == 'targets':
            __l_drug['targets'] = []

    elif checkParent('targets', stack):
        if tag == 'target':
            argDict = parseArgs(args)
            __l_target = {'partner':int(argDict['partner'])}
            try:
                __listappend(__targets[__l_target['partner']], __l_drug['drugbank-id'])
            except KeyError:
                __targets[__l_target['partner']] = [__l_drug['drugbank-id']]

    elif checkParent('target', stack):
        if tag == 'known-action':
            startCapture()

def handleEnd(tag, stack):
    global __tmp_drugs, __l_target, __excluded_targets

    if checkParent('drug', stack):
        if tag == 'drugbank-id':
            drug_bank_id = endCapture()
            __l_drug['drugbank-id'] = drug_bank_id
        elif tag == 'name':
            drug_name = endCapture()
            __l_drug['name'] = drug_name

    elif tag == 'drug' and checkParent('drugs', stack):
        __drugs[__l_drug['drugbank-id']] = __l_drug

    elif checkParent('targets', stack):
        if tag == 'target' and 'known-action' in __l_target.keys() and __l_target['known-action'] == 'yes':
            __listappend(__l_drug['targets'], __l_target)
        elif tag == 'target' and 'known-action' in __l_target.keys() and __l_target['known-action'] != 'yes':
            __excluded_targets += 1
    elif checkParent('target', stack):
        if tag == 'known-action':
            __l_target['known-action'] = endCapture()

def parseArgs(argString):
    args = {}
    if argString == None:
        return args
    
    items = __strsplit(__strstrip(argString), '"')
    
    next_type = 0
    lvar = ""
    for item in items:
        if next_type == 0:
            item = __strstrip(item)
            
            if item == '':
                continue
            if __strendswith(item, '='):
                next_type = 1
                item = item[:-1]
            item = __strstrip(item)

            args[item] = ""
            lvar = item
        elif next_type == 1:
            args[lvar] = item
            next_type = 0
    return args



def parseTag(tag):
    tag = __strstrip(tag)
    
    if __strstartswith(tag, "/"):
        return END, __strstrip(tag[1:]), None
    elif __strendswith(tag, "/"):
        tag = tag[:-1]
        name_end = __strfind(tag, " ")
        if name_end > -1:
            return CLOSURE, tag[:name_end], tag[name_end+1:]
        return CLOSURE, tag, None
    else:
        name_end = __strfind(tag, " ")
        if name_end > -1:
            return START, tag[:name_end], tag[name_end+1:]
        return START, tag, None

def parseLine(line, currentDrug, xmlStack):
    tag_iterator = re.finditer("<(.*?)>", line)

    last_index = 0
    
    next_element = tag_iterator.next
    
    try:
        while(True):
            pseudotag = next_element()
            ttype, tag, args = parseTag(pseudotag.group(1))
            
            tag_start = pseudotag.start(1)
            tag_end = pseudotag.end(1)
            
            capture(line[last_index:tag_start-1])
            last_index = tag_end + 1
            
            if ttype == CLOSURE:
                __listappend(xmlStack, tag)
                handleStart(tag, args, xmlStack)
                handleEnd(tag, xmlStack)
                __listpop(xmlStack)
            elif ttype == END:
                if xmlStack[-1] == tag:
                    handleEnd(tag, xmlStack)
                    __listpop(xmlStack)
                else:
                    raise SyntaxError("xml parsing error, incorrect nested tag: '%s'"
                            % (tag))
            else:
                __listappend(xmlStack, tag)
                handleStart(tag, args, xmlStack)

    except StopIteration:
        pass

class ProgressBar:
    def __init__(self):
        self.maximum = 100.0
        self.minimum = 0.0
        self.val = 0
        self.barwidth = 50

    def sysout(self):
        sys.stdout.write("\r")
        percentage = 100.0 * (self.val - self.minimum) / (self.valRange)
        
        num_bars = int(self.barwidth * percentage / 100.0)
        bars = '=' * num_bars
        if num_bars < self.barwidth:
            bars += ">"
        
        sys.stdout.write(("Progress: %6s |%-"+str(self.barwidth)+"s| ") % ("%.2f" % ( percentage ), bars))

    def setBarWidth(self, barwidth):
        self.barwidth = int(barwidth)

    def setProgress(self, progress):
        self.val = float(progress)

    def updateProgress(self, progress):
        self.val = float(progress)
        self.sysout()

    def setMininum(self, minval):
        self.minimum = float(minval)
        self.valRange = self.maximum - self.minimum

    def setMaximum(self, maxval):
        self.maximum = float(maxval)
        self.valRange = self.maximum - self.minimum

def saveCSV(filename):
    pass

def loadXML(filename):
    xmlfile = open(filename, 'r')
    
    stack = []
    
    cDrug = 0
    

    pbar = ProgressBar()
    pbar.setMaximum(float(2028522))
    lnum = 0
    
    for line in xmlfile:
        if lnum % 2500 == 0:
            pbar.updateProgress(lnum)
        lnum += 1
        cDrug = parseLine(line, cDrug, stack)
    
    pbar.updateProgress(lnum)
    print "Done"
    print "Targets excluded:           ", __excluded_targets
    print "Drugs loaded:               ", len(__drugs)
    print "Targets loaded:             ", len(__targets)

if __name__ == "__main__":
#    import psyco
#    psyco.full()
    geneDB.init(os.sep.join(["data","hgnc","hgnc_symbols.txt"]))
    
    print "Loading drug bank database from XML..."
    loadXML(os.sep.join(["data","drugbank","drugbank.xml"]))
    mapTargetNames(os.sep.join(["data","drugbank","target_links.csv"]))
    
    print __targets[1]
    print __targets[3]
    
    print "Named targets:   ", len(__target_names)
    print "Unnamed targets: ", __targets_unnamed

    print "Loading DrugBank drug catalogue..."
    drugCat.initDruglist(os.sep.join(["data","drugbank","drug_links.csv"]))
    print "Loading DrugBank drug target catalogue..."
    drugCat.initTargets(os.sep.join(["data","drugbank","target_links.csv"]),
            os.sep.join(["data","drugbank","all_target_protein.fasta"]),1,1)
    
    print "Drugcatalogue genes: ", len(drugCat.__geneSet)
    print "Overlaping genes:    ", len(__geneSet & drugCat.__geneSet)

# DB00191; DB00193; DB00226; DB00234; DB00285
    print __drugDict['slc6a2']
