import os
import re
import sys
__drugs = {}
__targets = {}
__excluded_targets = 0

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
            __l_target = {'partner':argDict['partner']}
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
    
    print "Loading drug bank database from XML..."

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
    print "Targets excluded: ", __excluded_targets
    print "Drugs loaded:     ", len(__drugs)
    print "Targets loaded:   ", len(__targets)

if __name__ == "__main__":
#    import psyco
#    psyco.full()

    loadXML(os.sep.join(["data","drugbank","drugbank.xml"]))
