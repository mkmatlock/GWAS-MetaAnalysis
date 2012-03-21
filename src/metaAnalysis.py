from geneList import *
import os,sys
import geneVerifier as geneDB
import gwasCatalog as gwasDB
import drugbankCatalog as drugDB
import math
import scipy.stats as stats
import htmltools
<<<<<<< HEAD
=======
from markup import oneliner
>>>>>>> html_reporting_refactor

__ENABLE_GENE_VERIFICATION = 0
__ENABLE_GENE_UPDATES = 0
__EXCLUDE_PSEUDOGENES = 0
__INCLUDE_MAPPED_GENES = 0
<<<<<<< HEAD

def getGeneListsByTrait(genes, pfilter):
=======
__EXCLUDE_OLFACTORY_PROTEINS = 1


__traitMetaAnalysis = {}

def computeTraitChiSquares(genes, pfilter):
>>>>>>> html_reporting_refactor
    
    traitSet = set([])
    for gene in genes:
        if gene in gwasDB.__traitDict:
            for trait in gwasDB.__traitDict[gene]:
                traitSet.add(trait)
    
    traitChi = {}
    
    for trait in traitSet:
        traitGenes = gwasDB.getGenesForTrait(trait,pfilter)
        
        listA = traitGenes & genes
        listC = traitGenes - genes
        
        a = len(listA)
        b = len(genes - traitGenes)
        c = len(listC)
        d = len(geneDB.__approved_symbols - (traitGenes | genes))
        
<<<<<<< HEAD
        chisq = geneUtils.contingentChiSquare(a,b,c,d)
        pvalue = stats.chisqprob(chisq, 1)
        
        traitChi[trait] = (a, chisq, pvalue, len(traitGenes))
    
    return traitChi

def writeGenePage(filename, title, desc, commonGenes):
    genepage = htmltools.createPage(title)
    
    total, geneTable = getGWASFrequencyTable(commonGenes)
=======
        chisq, oddsratio, kappa = geneUtils.contingentChiSquare(a,b,c,d)
        pvalue = stats.chisqprob(chisq, 1)
        
        traitChi[trait] = (a, chisq, oddsratio, kappa, pvalue, len(traitGenes))
    
    return traitChi
    

def writeGenePage(output_dir, filename, title, desc, total, geneTable):
    genepage = htmltools.createPage(title)

>>>>>>> html_reporting_refactor
    htmltools.pageDescription(genepage, desc)
    
    newTable=[]
    for row in geneTable:
<<<<<<< HEAD
        if os.path.exists("results\\html\\genelists"+os.sep+row[0]+".html"):
            newTable.append(["<a href=\"genelists/%s.html\">%s</a>" % (row[0], geneDB.__original_names[row[0]]), row[1]])
        else:
            newTable.append([geneDB.__original_names[row[0]], row[1]])
    
    htmltools.createTable(genepage, newTable, ["Gene","#Associated Traits"], "genelisthead", None, ["genecol","traitcol"],"genetable",None)
    
    htmltools.endPage(genepage)
    htmltools.savePage(genepage, filename)

=======
        drugcount = 0
        if row[0] in drugDB.__drugDict:
            drugcount = len(drugDB.__drugDict[row[0]])

        if os.path.exists(os.sep.join([output_dir,"genelists",row[0]+".html"])):
            newTable.append(["<a href=\"genelists/%s.html\">%s</a>" % (row[0], geneDB.__original_names[row[0]]), row[1], drugcount])
        else:
            newTable.append([geneDB.__original_names[row[0]], row[1], drugcount])
    
    htmltools.createTable(genepage, newTable, ["Gene","#Associated Traits", "#Targeted Drugs"], "genelisthead", None, ["genecol","traitcol", "drugcol"],"genetable",None)
    
    htmltools.endPage(genepage)
    htmltools.savePage(genepage, filename)
 
>>>>>>> html_reporting_refactor
def writeTraitPage(filename, title, desc, commonGenes, pfilter_cutoff):
    traitpage = htmltools.createPage(title)
    htmltools.pageDescription(traitpage, desc)
    
<<<<<<< HEAD
    traitChi = getGeneListsByTrait(commonGenes, pfilter_cutoff)
=======
    traitChi = computeTraitChiSquares(commonGenes, pfilter_cutoff)
>>>>>>> html_reporting_refactor
    
    traitTable = []
    
    for trait in traitChi:
<<<<<<< HEAD
        (cnt,chisq,pvalue,numgenes)=traitChi[trait]
                
        translate = trait.replace(" ","_").replace("/", " or ").replace("\\", " or ")
        
        if len(trait) > 40:
            trait = trait[:40] + "..."
        alink = "<a href=\"traitlists/%s.html\">%s</a>" % (translate, trait)
        
        entry = [alink, cnt, numgenes, "%.2f" % (chisq), "%.7f" % (pvalue)]
=======
        (cnt, chisq, oddsratio, kappa, pvalue,numgenes)=traitChi[trait]
        
        translate = trait.replace(" ","_").replace("/", " or ").replace("\\", " or ")
        
        if len(trait) > 38:
            trait = trait[:35] + "..."
        alink = "<a href=\"traitlists/%s.html\">%s</a>" % (translate, trait)
        
        entry = [alink, cnt, numgenes, "%.2f" % (chisq), "%.7f" % (pvalue), "%.1f" % (oddsratio), "%.4f" % (kappa)]
>>>>>>> html_reporting_refactor
        traitTable.append(entry)
        
    traitTable = sorted(traitTable, key=lambda item: -item[1])
    
<<<<<<< HEAD
    htmltools.createTable(traitpage, traitTable, ["Disease/Trait", "# RE Genes", "# Trait Genes", "chi-square", "p-value"], "traitlisthead", None, ["traitcol", "recol", "genecol", "chicol","pcol"], "traittable", None)
    htmltools.endPage(traitpage)
    htmltools.savePage(traitpage, filename)


def printList(list, columns = 4, width=30):
    i = 0
    while i < len(list):
        
        arg_string = ""
        arg_list = []
        
        j = i
        while j < i+columns and j < len(list):
            arg_string += "%-" + str(width) + "s "
            
            item = list[j]
            
            if (len(item) > width-1):
                item = item[:width-5] + "..."
            
            arg_list.append(item)
            j+=1
        
        print arg_string % tuple(arg_list)
        
        i+=columns
    
=======
    htmltools.createTable(traitpage, traitTable, ["Disease/Trait", "# RE Genes", "# Trait Genes", "chi-square", "p-value", "odds ratio", "kappa"], 
            "traitlisthead", None, ["traitcol", "recol", "genecol", "chicol", "pcol", "oddscol", "kappacol"], "traittable", None)
    htmltools.endPage(traitpage)
    htmltools.savePage(traitpage, filename)
    
def computeTraitDrugLists(RE_genes, drug_genes, pfilter_cutoff):
    
    traitSet = set(gwasDB.__studyByTrait.keys())
    
    for trait in traitSet:
        traitGenes = gwasDB.getGenesForTrait(trait, pfilter_cutoff)
        
        if len(traitGenes) == 0:
            continue
        
        RE_drugs = []
        other_drugs = []
        drug_counts_by_gene = {}
        
        for gene in traitGenes & RE_genes:
            drug_counts_by_gene[gene] = 0
            if gene in drugDB.__drugDict:
                for drug in drugDB.__drugDict[gene]:
                    RE_drugs.append(drug)
                drug_counts_by_gene[gene] += len(drugDB.__drugDict[gene])

        for gene in ((traitGenes & drug_genes) - RE_genes):
            drug_counts_by_gene[gene] = 0
            if gene in drugDB.__drugDict:
                for drug in drugDB.__drugDict[gene]:
                    other_drugs.append(drug)
                drug_counts_by_gene[gene] += len(drugDB.__drugDict[gene])
        __traitMetaAnalysis[trait]['RE_drugs'] = RE_drugs
        __traitMetaAnalysis[trait]['other_drugs'] = other_drugs
        __traitMetaAnalysis[trait]['drug_counts'] = drug_counts_by_gene


def computeTraitGeneLists(RE_genes, drug_genes, pfilter_cutoff):
    
    traitSet = set(gwasDB.__studyByTrait.keys())
    
    for trait in traitSet:
        traitGenes = gwasDB.getGenesForTrait(trait, pfilter_cutoff)
        
        if len(traitGenes) == 0: 
            continue
        
        __traitMetaAnalysis[trait] = {}
        
        
        RE = []
        for gene in traitGenes & RE_genes:
            
            count = 0
            for studyId in gwasDB.__studyByTrait[trait]:
                if studyId in gwasDB.__studyGenes:
                    if gene in gwasDB.__studyGenes[studyId]:
                        count+=1
            RE.append((gene, count))
        
        __traitMetaAnalysis[trait]['RE'] = RE
        
        
        
        drug = []
        for gene in traitGenes & drug_genes:
            
            count = 0
            for studyId in gwasDB.__studyByTrait[trait]:
                if studyId in gwasDB.__studyGenes:
                    if gene in gwasDB.__studyGenes[studyId]:
                        count+=1
            drug.append((gene, count))
        
        __traitMetaAnalysis[trait]['drugbank'] = drug
        
        
            
        other = []
        for gene in traitGenes - RE_genes - drug_genes:
            
            count = 0
            for studyId in gwasDB.__studyByTrait[trait]:
                if studyId in gwasDB.__studyGenes:
                    if gene in gwasDB.__studyGenes[studyId]:
                        count+=1
            other.append((gene,count))
        
        __traitMetaAnalysis[trait]['other'] = other
        
        
        
        a = len(traitGenes & RE_genes)
        b = len(RE_genes - traitGenes)
        c = len(traitGenes - RE_genes)
        d = len(geneDB.__approved_symbols - (traitGenes | RE_genes))
        
        chisq, oddsratio, kappa = geneUtils.contingentChiSquare(a,b,c,d)
        pvalue = stats.chisqprob(chisq, 1)
        
        __traitMetaAnalysis[trait]['RE_chi'] = (a, b, c, d, chisq, pvalue,
                oddsratio, kappa)
        
        a = len(traitGenes & drug_genes)
        b = len(drug_genes - traitGenes)
        c = len(traitGenes - drug_genes)
        d = len(geneDB.__approved_symbols - (traitGenes | drug_genes))
        
        chisq, oddsratio, kappa = geneUtils.contingentChiSquare(a,b,c,d)
        pvalue = stats.chisqprob(chisq, 1)
        
        __traitMetaAnalysis[trait]['drugbank_chi'] = (a, b, c, d, chisq,
                pvalue, oddsratio, kappa)
        
        __traitMetaAnalysis[trait]['geneset_size'] = len(traitGenes)
>>>>>>> html_reporting_refactor

def getGWASFrequencyTable(genes):
    output_list = []
    
    total = 0
    for geneSym in genes:
        if geneSym in gwasDB.__traitDict:
            l = len(gwasDB.__traitDict[geneSym])
        else:
            l = 0
        total += l
        output_list.append((geneSym, l))
        
    output_list = sorted(output_list, key=lambda item: -item[1])
    return total, output_list
<<<<<<< HEAD
    
    
def writeGeneListToFileWithGWASFreq(genes,filename):
    total, output_list = getGWASFrequencyTable(genes)
    
    ofile = open(filename,'w')
    ofile.write("%-20s  %d\n" % ("Total",total))
    ofile.write("-------------------------\n")
    for item in output_list:
        ofile.write("%-20s  %d\n" % item)
    ofile.close()
    
def writeGeneListToFile(genes,filename):

    geneList = sorted(list(genes))
    ofile = open(filename,'w')
    for gene in geneList:
        ofile.write(geneDB.__original_names[gene]+"\n")
    ofile.close()

def writeTraitFrequenciesToFile(genes, filename, geneListDir, pfilter):
        
    ofile = open(filename,'w')
    traitSet = set([])
    for gene in genes:
        if gene in gwasDB.__traitDict:
            for trait in gwasDB.__traitDict[gene]:
                traitSet.add(trait)
    
    probabilities = []
    
    for trait in traitSet:
        traitGenes = gwasDB.getGenesForTrait(trait,pfilter)
        
        listA = traitGenes & genes
        listC = traitGenes - genes
        if geneListDir != 0:
            geneListFilename = geneListDir + os.sep + trait.replace(" ","_").replace("/", " or ").replace("\\", " or ") + ".txt"
            
            geneFile = open(geneListFilename,'w')
            
            geneFile.write("%-20s%d\n" % ("Rapidly Evolving:",len(listA)))
            geneFile.write("---------------------\n")
            
            slist = []
            for gene in listA:
                
                count = 0
                for studyId in gwasDB.__studyByTrait[trait]:
                    if studyId in gwasDB.__studyGenes:
                        if gene in gwasDB.__studyGenes[studyId]:
                            count+=1
                slist.append((gene,count))
                
            for item in sorted(slist, key=lambda item: -item[1]):
                geneFile.write("%-20s%d\n" % item)
            
            geneFile.write("\n\n")
            geneFile.write("%-20s%d\n" % ("Other Genes:",len(listC)))
            geneFile.write("---------------------\n")
            
            slist = []
            for gene in listC:
                
                count = 0
                for studyId in gwasDB.__studyByTrait[trait]:
                    if studyId in gwasDB.__studyGenes:
                        if gene in gwasDB.__studyGenes[studyId]:
                            count+=1
                slist.append((gene,count))
                
            for item in sorted(slist, key=lambda item: -item[1]):
                geneFile.write("%-20s%d\n" % item)
                
            geneFile.close()
            
            
        a = len(listA)
        b = len(genes - traitGenes)
        c = len(listC)
        d = len(geneDB.__approved_symbols - (traitGenes | genes))
        
        
        
        chisq = geneUtils.contingentChiSquare(a,b,c,d)
        pvalue = stats.chisqprob(chisq, 1)
        
        if len(trait) > 48:
            trait = trait[0:45] + "..."
        
        probabilities.append((trait, str(a), str(len(traitGenes)), "%.2f" % (chisq), "%.10f" % (pvalue)))
        
        #       tG    ntG
        #  RE    a      b
        # nRE    c      d
        
    
    probabilities = sorted(probabilities, key=lambda item: -int(item[1]))
    
    ofile.write("%-50s%5s%8s%15s%15s\n" % ("Disease/Trait", "#RE", "#Genes", "chi-square", "p-value"))
    ofile.write("---------------------------------------------------------------------------------------------\n")
    for entry in probabilities:
        ofile.write("%-50s%5s / %5s%15s%15s\n" % entry)
        
    ofile.close()
    return traitSet
=======

    
def createGeneListTable(page, genes, verify=None):
    output_list = []
    
    total = 0
    for pair in genes:
        geneSym, l = pair
        total += l
        if geneSym in drugDB.__drugDict:
            drugs = len(drugDB.__drugDict[geneSym])
        else:
            drugs = 0
        if verify==None or geneSym in verify:
            output_list.append(["<a href=\"../genelists/%s.html\" > %s </a>" % (geneSym, geneDB.__original_names[geneSym]), l, drugs])
        else:
            output_list.append([geneDB.__original_names[geneSym], l, drugs])
        
    output_list = sorted(output_list, key=lambda item: -item[1])
    htmltools.createTable(page, output_list, ["Gene", "#Associated Traits", "#Targeted Drugs"], "genelisthead", None, ["genecol","traitcol", "drugcol"],"genetable",None)
    
    
    
def createTraitListingsHTML(traitListDir):
    
    traitSet = set(gwasDB.__studyByTrait.keys())
    
    for trait in traitSet:
        if trait not in __traitMetaAnalysis:
            continue
        
        traitMetadata = __traitMetaAnalysis[trait]
        RE_proteins = traitMetadata['RE']
        drug_proteins = traitMetadata['drugbank']
        other_proteins = traitMetadata['other']
        
        chi_RE = traitMetadata['RE_chi']
        chi_Drugbank = traitMetadata['drugbank_chi']
        
        traitListFilename = os.sep.join([traitListDir, trait.replace(" ","_").replace("/", " or ").replace("\\", " or ") + ".html"])
        
        traitpage = htmltools.createPage("Trait Summary: " + trait, css_file='../genereport.css')
        
        htmltools.pageDescription(traitpage, "Gene list overlap summary for trait: %s" % (trait))
        
        # two of these
        
        htmltools.createChiTable(traitpage, "Overlap with rapidly evolving genes:", "RE", "trait", chi_RE[0], chi_RE[1], chi_RE[2], 
                chi_RE[3], chi_RE[4], chi_RE[5], chi_RE[6], chi_RE[7])
        
        htmltools.createChiTable(traitpage, "Overlap with drugbank genes:",
                "Drugbank", "trait", chi_Drugbank[0], chi_Drugbank[1],
                chi_Drugbank[2], chi_Drugbank[3], chi_Drugbank[4],
                chi_Drugbank[5], chi_Drugbank[6], chi_Drugbank[7])
        
        traitpage.table.open(class_="invisible")
        traitpage.tr.open()

        traitpage.td.open()
        traitpage.div("Gene Lists:", class_="header")
        
        traitpage.div("Trait genes indicated as rapidly evolving: ", class_="description")
        createGeneListTable(traitpage, RE_proteins)
        
        traitpage.div("Trait genes associated with Drugbank targets: ", class_="description")
        createGeneListTable(traitpage, drug_proteins)
        
        traitpage.div("Other trait genes: ", class_="description")
        createGeneListTable(traitpage, other_proteins)
        
        traitpage.td.close()
        traitpage.td.open()

        traitpage.div("Drug Lists:", class_="header")
        
        druglistlen = len(__traitMetaAnalysis[trait]['RE_drugs'])
        traitpage.div("%d drugs targeting associated rapidly evolving proteins" % (druglistlen), class_="description")
        
        traitpage.div.open(class_="druglist")
        traitpage.ul.open()
        for drug in __traitMetaAnalysis[trait]['RE_drugs']:
            link = "http://www.drugbank.ca/drugs/%s" % (drug)
            if drug not in drugDB.__drugs:
                traitpage.li(oneliner.a(drug, href=link))
            elif drugDB.__drugs[drug][1] != None:
                traitpage.li("<a href = \"%s\">%s</a> [<a href=\"%s\">alt</a>]" % (link, drugDB.__drugs[drug][0], drugDB.__drugs[drug][1]))
            else:
                traitpage.li(oneliner.a(drugDB.__drugs[drug][0], href=link))
        traitpage.ul.close()
        traitpage.div.close()
        
        druglistlen = len(__traitMetaAnalysis[trait]['other_drugs'])
        traitpage.div("%d Drugs targeting other proteins" % (druglistlen), class_="description")
        

        traitpage.div.open(class_="druglist")
        traitpage.ul.open()
        for drug in __traitMetaAnalysis[trait]['other_drugs']:
            link = "http://www.drugbank.ca/drugs/%s" % (drug)
            if drug not in drugDB.__drugs:
                traitpage.li(oneliner.a(drug, href=link))
            elif drugDB.__drugs[drug][1] != None:
                traitpage.li("<a href = \"%s\">%s</a> [<a href=\"%s\">alt</a>]" % (link, drugDB.__drugs[drug][0], drugDB.__drugs[drug][1]))
            else:
                traitpage.li(oneliner.a(drugDB.__drugs[drug][0], href=link))
        traitpage.ul.close()
        traitpage.div.close()


        traitpage.td.close()
        traitpage.tr.close()
        traitpage.table.close()
        
        htmltools.savePage(traitpage, traitListFilename)
    
def createGeneListingsHTML(geneListDir):
    for gene in gwasDB.__geneSet:
        
        genePage = htmltools.createPage("Gene Summary: " + geneDB.__original_names[gene], css_file='../genereport.css')
        
        
        # Create the disease trait tables
        traits = []
        for trait in gwasDB.__traitDict[gene]:
            traits.append(trait)
            
        traits = sorted(traits, key=lambda trait: -__traitMetaAnalysis[trait]['RE_chi'][0])
        
        traitTable = []
        for trait in traits:
            cnt       = len(__traitMetaAnalysis[trait]['RE'])
            chisq     = __traitMetaAnalysis[trait]['RE_chi'][4]
            pvalue    = __traitMetaAnalysis[trait]['RE_chi'][5]
            numgenes  = __traitMetaAnalysis[trait]['geneset_size']
            translate = trait.replace(" ","_").replace("/", " or ").replace("\\", " or ")
            
            if len(trait) > 38:
                trait = trait[:35] + "..."
            traitTable.append(["<a href=\"../traitlists/%s.html\">%s</a>" % (translate,trait), cnt, numgenes, "%.2f" % (chisq), "%.7f" % (pvalue)])
            
        genePage.div("Gene %s, total traits: %d" % (geneDB.__original_names[gene],len(traitTable)), class_="header")
        
        htmltools.createTable(genePage, traitTable, ["Disease/Trait", "#RE Genes", "#Trait Genes", "Chi-square", "p-value"], "traitlisthead", None, ["traitcol","recol", "genecol", "chicol","pcol"], "traittable", None)
        
        # Create drug bank links
        
        
        if gene not in drugDB.__drugDict:
            genePage.div("No drugs target gene %s" % (geneDB.__original_names[gene]), class_="header")
        else:
            drugbank_size = len(drugDB.__drugDict[gene])
            
            genePage.div("%d drugs targeting gene %s" % (drugbank_size, geneDB.__original_names[gene]), class_="header")
            
            genePage.div.open(class_="druglist")
            genePage.ul.open()
            for drug in drugDB.__drugDict[gene]:
                link = "http://www.drugbank.ca/drugs/%s" % (drug)
                if drug not in drugDB.__drugs:
                    genePage.li(oneliner.a(drug, href=link))
                elif drugDB.__drugs[drug][1] != None:
                    genePage.li("<a href = \"%s\">%s</a> [<a href=\"%s\">alt</a>]" % (link, drugDB.__drugs[drug][0], drugDB.__drugs[drug][1]))
                else:
                    genePage.li(oneliner.a(drugDB.__drugs[drug][0], href=link))
            genePage.ul.close()
            genePage.div.close()
            
        
        htmltools.savePage(genePage, os.sep.join([geneListDir, gene + ".html"]))
>>>>>>> html_reporting_refactor

if __name__ == "__main__":
    
    pfilter_cutoff = 0.05
    trait_exclude_file = 0
    gene_frequency_filter=1
    skip_listings=0
<<<<<<< HEAD
=======
    output_dir = os.sep.join(["results", "html"])
    included_studies = [1,2,3,4,5]
>>>>>>> html_reporting_refactor
    for i in xrange(1, len(sys.argv)):
        if sys.argv[i] == "--pfilter":
            pfilter_cutoff = float(sys.argv[i+1])
        elif sys.argv[i] == "--verify":
            __ENABLE_GENE_VERIFICATION = 1
        elif sys.argv[i] == "--update":
            __ENABLE_GENE_UPDATES = 1
        elif sys.argv[i] == "--inc-map":
            __INCLUDE_MAPPED_GENES = 1
        elif sys.argv[i] == "--exclude-traits":
            trait_exclude_file = sys.argv[i+1]
        elif sys.argv[i] == "--freq":
            gene_frequency_filter = int(sys.argv[i+1])
        elif sys.argv[i] == "--skip-lists":
            skip_listings=1
        elif sys.argv[i] == "--rm-pseudo":
            __EXCLUDE_PSEUDOGENES = 1
<<<<<<< HEAD
=======
        elif sys.argv[i] == "--rm-study":
            for item in sys.argv[i+1].split(","):
                item = item.strip()
                included_studies.remove(int(item))
        elif sys.argv[i] == "--output-dir":
            output_dir = os.sep.join(["results", sys.argv[i+1]])
>>>>>>> html_reporting_refactor
        elif sys.argv[i].startswith("--"):
            print "Invalid arg:", sys.argv[i]
            sys.exit(1)
    
<<<<<<< HEAD
    # clean the geneList dir
    print "Cleaning gene lists..."
    geneListFolder = "results"+os.sep+"geneLists"
    for filename in os.listdir(geneListFolder):
        file_path = os.path.join(geneListFolder, filename)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except:
            pass
    
    geneDB.init("data\\genelists\\hgnc_symbols.txt",__EXCLUDE_PSEUDOGENES)
    
    print "\nLoading rapidly evolving proteins..."
    studyGenes = loadEvolutionaryGenes("data\\genelists\\positive_selection_gene_lists.txt",__ENABLE_GENE_VERIFICATION,__ENABLE_GENE_UPDATES,gene_frequency_filter)
    
    print "\nLoading GWAS catalogue..."
    gwasDB.init("data\\gwas\\gwascatalog.txt",__ENABLE_GENE_VERIFICATION, __ENABLE_GENE_UPDATES, __INCLUDE_MAPPED_GENES,trait_exclude_file,pfilter_cutoff)
    
    print "\nLoading DrugBank drug catalogue..."
    drugDB.initTargets("data\\drugbank\\drug_links.csv")
    print "\nLoading DrugBank drug target catalogue..."
    drugDB.initTargets("data\\drugbank\\target_links.csv", "data\\drugbank\\all_target_protein.fasta",__ENABLE_GENE_VERIFICATION,__ENABLE_GENE_UPDATES)
=======
    print "\nInitializing HGNC Database..."
    geneDB.init(os.sep.join(["data","hgnc","hgnc_symbols.txt"]),__EXCLUDE_PSEUDOGENES)
    
    print "\nLoading rapidly evolving proteins..."
    studyGenes = loadEvolutionaryGenes(os.sep.join(["data","genelists","positive_selection_gene_lists.txt"]),__ENABLE_GENE_VERIFICATION,__ENABLE_GENE_UPDATES,gene_frequency_filter, included_studies)
    
    if __EXCLUDE_OLFACTORY_PROTEINS:
        print "\nExcluding olfactory proteins..."
        olfactoryProteins = pyCSV()
        olfactoryProteins.load(os.sep.join(["data","genelists","olfactory_genes.csv"]),"\t")
        
        olfactoryGenes = set([item.lower() for item in geneUtils.columnToList(olfactoryProteins, 1, 1)])
        
        for gene in olfactoryGenes:
            if gene in studyGenes:
                studyGenes.remove(gene)
    
    print "\nLoading GWAS catalogue..."
    gwasDB.init(os.sep.join(["data","gwas","gwascatalog.txt"]),__ENABLE_GENE_VERIFICATION, __ENABLE_GENE_UPDATES, __INCLUDE_MAPPED_GENES,trait_exclude_file,pfilter_cutoff)
    print "\nLoading DrugBank drug catalogue..."
    drugDB.initDruglist(os.sep.join(["data","drugbank","drug_links.csv"]))
    print "\nLoading DrugBank drug target catalogue..."
    drugDB.initTargets(os.sep.join(["data","drugbank","target_links.csv"]), os.sep.join(["data","drugbank","all_target_protein.fasta"]),__ENABLE_GENE_VERIFICATION,__ENABLE_GENE_UPDATES)
>>>>>>> html_reporting_refactor
    
    commonGenes = studyGenes & gwasDB.__geneSet
    cgWithoutGWAS = studyGenes - gwasDB.__geneSet
    gwasWithoutCG = gwasDB.__geneSet - studyGenes
    
    
<<<<<<< HEAD
    if not os.path.exists("results" + os.sep + "html"):
        os.mkdir("results" + os.sep + "html")
        os.mkdir("results" + os.sep + "html" + os.sep + "genelists")
        os.mkdir("results" + os.sep + "html" + os.sep + "traitlists")
    
    # make index.html
    
=======
    # chi-square contingency table for gwas and RE
    
    a1 = len(commonGenes)
    b1 = len(gwasWithoutCG)
    c1 = len(cgWithoutGWAS)
    d1 = len( geneDB.__approved_symbols - (studyGenes | gwasDB.__geneSet) )
    chisq1, oddsratio1, kappa1 = geneUtils.contingentChiSquare(a1,b1,c1,d1)
    pvalue1 = stats.chisqprob(chisq1, 1)
    
    
    
    commonDrugTargets = drugDB.__geneSet & studyGenes
    
    a2 = len( commonDrugTargets )
    b2 = len( drugDB.__geneSet - studyGenes )
    c2 = len( studyGenes - drugDB.__geneSet )
    d2 = len( geneDB.__approved_symbols - ( studyGenes | drugDB.__geneSet ) )
    chisq2, oddsratio2, kappa2  = geneUtils.contingentChiSquare(a2,b2,c2,d2)
    pvalue2 = stats.chisqprob(chisq2, 1)
    
    
    
    
    gwas_drugbank_overlap = gwasDB.__geneSet & drugDB.__geneSet
    
    a3 = len( gwas_drugbank_overlap )
    b3 = len( drugDB.__geneSet - gwasDB.__geneSet )
    c3 = len( gwasDB.__geneSet - drugDB.__geneSet )
    d3 = len( geneDB.__approved_symbols - ( drugDB.__geneSet | gwasDB.__geneSet ) )
    chisq3, oddsratio3, kappa3  = geneUtils.contingentChiSquare(a3,b3,c3,d3)
    pvalue3 = stats.chisqprob(chisq3, 1)
    
    
    # start preparing HTML reports
    
    print "\n------------------------------\nPreparing HTML Reports..."
    
    if not os.path.exists("results"):
        os.mkdir("results")
    if not os.path.exists(os.sep.join(["results","log"])):
        os.mkdir(os.sep.join(["results", "log"]))
    if not os.path.exists(os.sep.join([output_dir])):
        os.mkdir(os.sep.join([output_dir]))
        os.mkdir(os.sep.join([output_dir, "genelists"]))
        os.mkdir(os.sep.join([output_dir, "traitlists"]))
        os.mkdir(os.sep.join([output_dir, "DAVID"]))
        
    print "\nCopying support files..."
        
    os.system("copy %s %s" % (os.sep.join(["html", "genereport.css"]), output_dir ))
    os.system("copy %s %s" % (os.sep.join(["html", "DAVID", str(gene_frequency_filter), "*"]), os.sep.join([output_dir, "DAVID"])))
    
    # make index.html
    
    print "Making index.html..."
    
    # Report data set properties
    
>>>>>>> html_reporting_refactor
    indexpage = htmltools.createPage("Rapidly Evolving Gene Meta-Analysis Report")
    
    indexpage.div("Summary of loaded information:", class_="header")
    
    indexpage.div("Rapidly Evolving genes loaded:                  %d".replace(" ", "&nbsp;") % (len(studyGenes)),class_="console")
    indexpage.div("Drugbank drug target proteins loaded:           %d".replace(" ", "&nbsp;") % (len(drugDB.__geneSet)),class_="console")
    indexpage.div("GWAS Genes indicated in disease studies:        %d".replace(" ", "&nbsp;") % (len(gwasDB.__geneSet)),class_="console")
    
<<<<<<< HEAD
    a = len(commonGenes)
    b = len(gwasWithoutCG)
    c = len(cgWithoutGWAS)
    d = len( geneDB.__approved_symbols - (studyGenes | gwasDB.__geneSet) )
    chisq = geneUtils.contingentChiSquare(a,b,c,d)
    pvalue = stats.chisqprob(chisq, 1)
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "GWAS Overlap with Rapidly Evolving Geneset:", "GWAS", "RE", a, b, c, d, chisq, pvalue)
=======
    
    # report gwas, RE overlap chi matrix
    
    htmltools.createChiTable(indexpage, "GWAS Overlap with Rapidly Evolving Geneset:", "GWAS", "RE", a1, b1, c1, d1, chisq1, pvalue1,
            oddsratio1, kappa1)
>>>>>>> html_reporting_refactor
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="gwas_traits.html")
    indexpage.a("Gene Listing", href="gwas_genes.html")
    
<<<<<<< HEAD
    if gene_frequency_filter == 1:
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_gwas_common.xhtml")
        
        
    indexpage.div.close()
    
    indexpage.div.close()
    
    print "\n\n----------GWAS Results----------"
    print "Test Matrix:       %5s   %5s" % ("RE", "nRE")
    print "             GWAS  %5s   %5s" % (str(a), str(b))
    print "            nGWAS  %5s   %5s" % (str(c), str(d))
    print "Chi statistic:    %.2f" % (chisq)
    print "P-value:          %f" % (pvalue)
    print "--------------------------------"
    
    
    a = len( studyGenes & drugDB.__geneSet )
    b = len( drugDB.__geneSet - studyGenes )
    c = len( studyGenes - drugDB.__geneSet )
    d = len( geneDB.__approved_symbols - ( studyGenes | drugDB.__geneSet ) )
    chisq  = geneUtils.contingentChiSquare(a,b,c,d)
    pvalue = stats.chisqprob(chisq, 1)
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "Drugbank Overlap with Rapidly Evolving Geneset:", "Drug", "RE", a, b, c, d, chisq, pvalue)
=======
    if os.path.exists(os.sep.join([output_dir, "DAVID", "david_gwas_common.xhtml"])):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_gwas_common.xhtml")
        
    indexpage.div.close()
    
    
    # report drugbank, RE overlap chi matrix
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "Drugbank Overlap with Rapidly Evolving Geneset:", "Drugbank", "RE", a2, b2, c2, d2, chisq2, pvalue2, oddsratio2, kappa2)
>>>>>>> html_reporting_refactor
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="overlap_traits.html")
    indexpage.a("Gene Listing", href="drugbank_genes.html")
    
<<<<<<< HEAD
    if gene_frequency_filter == 1:
=======
    if os.path.exists(os.sep.join([output_dir, "DAVID", "david_drugbank_common.xhtml"])):
>>>>>>> html_reporting_refactor
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_drugbank_common.xhtml")
    
    indexpage.div.close()
    
<<<<<<< HEAD
    indexpage.div.close()
    
    print "\n\n--------Drugbank Results--------"
    print "Test Matrix:       %5s   %5s" % ("RE", "nRE")
    print "             Drug  %5s   %5s" % (str(a), str(b))
    print "            nDrug  %5s   %5s" % (str(c), str(d))
    print "Chi statistic:    %.2f" % (chisq)
    print "P-value:          %f" % (pvalue)
    print "--------------------------------"
    
    commonDrugTargets = drugDB.__geneSet & studyGenes
    
    gwas_drugbank_overlap=gwasDB.__geneSet & drugDB.__geneSet
    a = len( gwas_drugbank_overlap )
    b = len( drugDB.__geneSet - gwasDB.__geneSet )
    c = len( gwasDB.__geneSet - drugDB.__geneSet )
    d = len( geneDB.__approved_symbols - ( drugDB.__geneSet | gwasDB.__geneSet ) )
    chisq  = geneUtils.contingentChiSquare(a,b,c,d)
    pvalue = stats.chisqprob(chisq, 1)
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "GWAS Overlap with Drugbank:", "Drugbank", "GWAS", a, b, c, d, chisq, pvalue)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="drugbank_traits.html", style="align:left")
    indexpage.a("Gene Listing", href="drugbank_gwas_genes.html", style="align:right")
    
    if gene_frequency_filter == 1:
=======
    
    # report drugbank, GWAS overlap chi matrix
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "GWAS Overlap with Drugbank:",
            "Drugbank", "GWAS", a3, b3, c3, d3, chisq3, pvalue3, 
            oddsratio3, kappa3)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="drugbank_traits.html")
    indexpage.a("Gene Listing", href="drugbank_gwas_genes.html")
    
    if os.path.exists(os.sep.join([output_dir, "DAVID", "david_drugbank_gwas_common.xhtml"])):
>>>>>>> html_reporting_refactor
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_drugbank_gwas_common.xhtml")
    
    indexpage.div.close()
    
<<<<<<< HEAD
    indexpage.div.close()
    
    
    
    print "\n\n--------Combined Results--------"
    print "Overlap (GWAS / Drugbank):        "
    print "Test Matrix:       %5s   %5s" % ("GWAS", "nGWAS")
    print "             Drug  %5s   %5s" % (str(a), str(b))
    print "            nDrug  %5s   %5s" % (str(c), str(d))
    print "Chi statistic:    %.2f" % (chisq)
    print "P-value:          %f" % (pvalue)
    print "--------------------------------"
    print ""
=======
>>>>>>> html_reporting_refactor
    
    overlap = commonGenes & studyGenes & drugDB.__geneSet
    
    indexpage.div.open(class_="links")
    indexpage.p.open()
    indexpage.add("Total overlap of drugbank, GWAS, and Rapidly Evolving geneset: ")
    indexpage.a(str(len(overlap)), href="all_genes.html")
    
<<<<<<< HEAD
    if gene_frequency_filter == 1: 
=======
    if os.path.exists(os.sep.join([output_dir, "DAVID", "david_all.xhtml"])):
>>>>>>> html_reporting_refactor
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_all.xhtml")
    
    indexpage.p.close()
<<<<<<< HEAD
    
    
    
    indexpage.div.close()
    
    print "Overlap (Everything):      %5s" % (str(len(overlap)))
    print "--------------------------------"
    
    
    writeGeneListToFile(gwas_drugbank_overlap, "results\\common_genes_gwas_drugbank.txt")
    
    writeGeneListToFile(overlap, "results\\common_genes_all.txt")
    writeGeneListToFileWithGWASFreq(overlap, "results\\common_genes_all_freq.txt")
    
    print "\n\nWriting common GWAS genes..."
    writeGeneListToFile(commonGenes, "results\\common_genes_gwas.txt")
    writeGeneListToFileWithGWASFreq(commonGenes, "results\\common_genes_gwas_freq.txt")
    
    print "\n\nWriting common GWAS gene trait frequencies..."
    gwas_traits = writeTraitFrequenciesToFile(studyGenes, "results\\gwas_trait_file.txt", "results\\geneLists", pfilter_cutoff)
        
    print "\n\nWriting exceptional genes..."
    writeGeneListToFile(cgWithoutGWAS, "results\\excepted_genes_gwas.txt")
    
    print "\n\nWriting common Drugbank genes..."
    writeGeneListToFile(commonDrugTargets, "results\\common_genes_drugbank.txt")
    
    drug_traits = writeTraitFrequenciesToFile(gwasDB.__geneSet & drugDB.__geneSet, "results\\drugbank_trait_file.txt", 0, pfilter_cutoff)
    
    print "\n\nWriting common Drugbank gene trait frequencies..."
    all_traits = writeTraitFrequenciesToFile(commonDrugTargets, "results\\all_trait_file.txt", 0, pfilter_cutoff)
    
    print "\n\n-----------------------------"
    print "GWAS Traits:        ", len(gwas_traits)
    print "Target Traits:      ", len(drug_traits)
    print "All traits:         ", len(all_traits)
    print "-----------------------------\n"
    printList(sorted(list(all_traits)), columns=4, width=35)
    
    htmltools.endPage(indexpage)
    htmltools.savePage(indexpage, "results"+os.sep+"html"+os.sep+"index.html")
=======
    indexpage.div.close()
    
    htmltools.endPage(indexpage)
    htmltools.savePage(indexpage, os.sep.join([output_dir, "index.html"]))
>>>>>>> html_reporting_refactor
    
    print ""
    
    if not skip_listings:
        print "Creating trait listings..."
<<<<<<< HEAD
        traitChi = htmltools.createTraitListingsHTML("results"+os.sep+"html"+os.sep+"traitlists", commonGenes, pfilter_cutoff)
        print "Creating gene listings..."
        htmltools.createGeneListingsHTML("results"+os.sep+"html"+os.sep+"genelists", gwasDB.__geneSet, traitChi)
=======
        computeTraitGeneLists(studyGenes, drugDB.__geneSet, pfilter_cutoff)
        computeTraitDrugLists(studyGenes, drugDB.__geneSet, pfilter_cutoff)
        createTraitListingsHTML(os.sep.join([output_dir, "traitlists"]))
        print "Creating gene listings..."
        createGeneListingsHTML(os.sep.join([output_dir,"genelists"]))
>>>>>>> html_reporting_refactor
    
    print "Creating referenced HTML gene and trait list reports..."
    
    # write trait frequency tables
<<<<<<< HEAD
    writeTraitPage("results\\html\\gwas_traits.html", "GWAS Traits for Rapidly Evolving Genes", "Listing of disease traits associated with genes in the rapidly evolving geneset", commonGenes, pfilter_cutoff)
    writeTraitPage("results\\html\\drugbank_traits.html", "GWAS Traits for Drugbank Genes", "Listing of disease traits associated with genes in the drugbank targets database", gwasDB.__geneSet & drugDB.__geneSet, pfilter_cutoff)
    writeTraitPage("results\\html\\overlap_traits.html", "GWAS Traits for overlap of RE and Drugbank", "Listing of disease traits associated with genes in the rapidly evolving geneset and the drugbank targets database", overlap, pfilter_cutoff)
    
    # write gene tables
    
    writeGenePage("results\\html\\gwas_genes.html", "Rapidly Evolving Genes found in GWAS Disease Studies", "Listing of rapidly evolving genes found in GWAS.", commonGenes)
    writeGenePage("results\\html\\drugbank_genes.html", "Rapidly Evolving Genes found in Drugbank", "Listing of rapidly evolving genes found in Drugbank.", commonDrugTargets)
    writeGenePage("results\\html\\drugbank_gwas_genes.html", "Drugbank targets found in GWAS Disease Studies", "Listing of Drugbank targets found in GWAS Disease Studies.", gwas_drugbank_overlap)
    writeGenePage("results\\html\\all_genes.html", "Rapidly Evolving Genes found in GWAS Disease Studies and Drugbank", "Listing of rapidly evolving genes found in both GWAS and Drugbank.", overlap)
    
    
    
    
=======
    writeTraitPage(os.sep.join([output_dir,"gwas_traits.html"]), "GWAS Traits for Rapidly Evolving Genes", "Listing of disease traits associated with genes in the rapidly evolving geneset", commonGenes, pfilter_cutoff)
    writeTraitPage(os.sep.join([output_dir,"drugbank_traits.html"]), "GWAS Traits for Drugbank Genes", "Listing of disease traits associated with genes in the drugbank targets database", gwasDB.__geneSet & drugDB.__geneSet, pfilter_cutoff)
    writeTraitPage(os.sep.join([output_dir,"overlap_traits.html"]), "GWAS Traits for overlap of RE and Drugbank", "Listing of disease traits associated with genes in the rapidly evolving geneset and the drugbank targets database", overlap, pfilter_cutoff)
    
    # write gene tables
    
    total, geneTable = getGWASFrequencyTable(commonGenes)
    writeGenePage(output_dir, os.sep.join([output_dir,"gwas_genes.html"]), "Rapidly Evolving Genes found in GWAS Disease Studies", "Listing of rapidly evolving genes found in GWAS.", total, geneTable)
    
    total, geneTable = getGWASFrequencyTable(commonDrugTargets)
    writeGenePage(output_dir, os.sep.join([output_dir,"drugbank_genes.html"]), "Rapidly Evolving Genes found in Drugbank", "Listing of rapidly evolving genes found in Drugbank.", total, geneTable)
    
    total, geneTable = getGWASFrequencyTable(gwas_drugbank_overlap)
    writeGenePage(output_dir, os.sep.join([output_dir,"drugbank_gwas_genes.html"]), "Drugbank targets found in GWAS Disease Studies", "Listing of Drugbank targets found in GWAS Disease Studies.", total, geneTable)
    
    total, geneTable = getGWASFrequencyTable(overlap)
    writeGenePage(output_dir, os.sep.join([output_dir,"all_genes.html"]), "Rapidly Evolving Genes found in GWAS Disease Studies and Drugbank", "Listing of rapidly evolving genes found in both GWAS and Drugbank.", total, geneTable)
    
    
    
    
>>>>>>> html_reporting_refactor
