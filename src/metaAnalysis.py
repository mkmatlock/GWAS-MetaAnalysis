from geneList import *
import os,sys
import geneVerifier as geneDB
import gwasCatalog as gwasDB
import drugbankCatalog as drugDB
import math
import scipy.stats as stats
import htmltools

__ENABLE_GENE_VERIFICATION = 0
__ENABLE_GENE_UPDATES = 0
__EXCLUDE_PSEUDOGENES = 0
__INCLUDE_MAPPED_GENES = 0

def getGeneListsByTrait(genes, pfilter):
    
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
        
        chisq = geneUtils.contingentChiSquare(a,b,c,d)
        pvalue = stats.chisqprob(chisq, 1)
        
        traitChi[trait] = (a, chisq, pvalue, len(traitGenes))
    
    return traitChi

def writeGenePage(filename, title, desc, commonGenes):
    genepage = htmltools.createPage(title)
    
    total, geneTable = getGWASFrequencyTable(commonGenes)
    htmltools.pageDescription(genepage, desc)
    
    newTable=[]
    for row in geneTable:
        if os.path.exists("results\\html\\genelists"+os.sep+row[0]+".html"):
            newTable.append(["<a href=\"genelists/%s.html\">%s</a>" % (row[0], geneDB.__original_names[row[0]]), row[1]])
        else:
            newTable.append([geneDB.__original_names[row[0]], row[1]])
    
    htmltools.createTable(genepage, newTable, ["Gene","#Associated Traits"], "genelisthead", None, ["genecol","traitcol"],"genetable",None)
    
    htmltools.endPage(genepage)
    htmltools.savePage(genepage, filename)

def writeTraitPage(filename, title, desc, commonGenes, pfilter_cutoff):
    traitpage = htmltools.createPage(title)
    htmltools.pageDescription(traitpage, desc)
    
    traitChi = getGeneListsByTrait(commonGenes, pfilter_cutoff)
    
    traitTable = []
    
    for trait in traitChi:
        (cnt,chisq,pvalue,numgenes)=traitChi[trait]
                
        translate = trait.replace(" ","_").replace("/", " or ").replace("\\", " or ")
        
        if len(trait) > 40:
            trait = trait[:40] + "..."
        alink = "<a href=\"traitlists/%s.html\">%s</a>" % (translate, trait)
        
        entry = [alink, cnt, numgenes, "%.2f" % (chisq), "%.7f" % (pvalue)]
        traitTable.append(entry)
        
    traitTable = sorted(traitTable, key=lambda item: -item[1])
    
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

if __name__ == "__main__":
    
    pfilter_cutoff = 0.05
    trait_exclude_file = 0
    gene_frequency_filter=1
    skip_listings=0
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
        elif sys.argv[i].startswith("--"):
            print "Invalid arg:", sys.argv[i]
            sys.exit(1)
    
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
    
    commonGenes = studyGenes & gwasDB.__geneSet
    cgWithoutGWAS = studyGenes - gwasDB.__geneSet
    gwasWithoutCG = gwasDB.__geneSet - studyGenes
    
    
    if not os.path.exists("results" + os.sep + "html"):
        os.mkdir("results" + os.sep + "html")
        os.mkdir("results" + os.sep + "html" + os.sep + "genelists")
        os.mkdir("results" + os.sep + "html" + os.sep + "traitlists")
    
    # make index.html
    
    indexpage = htmltools.createPage("Rapidly Evolving Gene Meta-Analysis Report")
    
    indexpage.div("Summary of loaded information:", class_="header")
    
    indexpage.div("Rapidly Evolving genes loaded:                  %d".replace(" ", "&nbsp;") % (len(studyGenes)),class_="console")
    indexpage.div("Drugbank drug target proteins loaded:           %d".replace(" ", "&nbsp;") % (len(drugDB.__geneSet)),class_="console")
    indexpage.div("GWAS Genes indicated in disease studies:        %d".replace(" ", "&nbsp;") % (len(gwasDB.__geneSet)),class_="console")
    
    a = len(commonGenes)
    b = len(gwasWithoutCG)
    c = len(cgWithoutGWAS)
    d = len( geneDB.__approved_symbols - (studyGenes | gwasDB.__geneSet) )
    chisq = geneUtils.contingentChiSquare(a,b,c,d)
    pvalue = stats.chisqprob(chisq, 1)
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "GWAS Overlap with Rapidly Evolving Geneset:", "GWAS", "RE", a, b, c, d, chisq, pvalue)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="gwas_traits.html")
    indexpage.a("Gene Listing", href="gwas_genes.html")
    
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
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="overlap_traits.html")
    indexpage.a("Gene Listing", href="drugbank_genes.html")
    
    if gene_frequency_filter == 1:
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_drugbank_common.xhtml")
    
    indexpage.div.close()
    
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
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_drugbank_gwas_common.xhtml")
    
    indexpage.div.close()
    
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
    
    overlap = commonGenes & studyGenes & drugDB.__geneSet
    
    indexpage.div.open(class_="links")
    indexpage.p.open()
    indexpage.add("Total overlap of drugbank, GWAS, and Rapidly Evolving geneset: ")
    indexpage.a(str(len(overlap)), href="all_genes.html")
    
    if gene_frequency_filter == 1: 
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_all.xhtml")
    
    indexpage.p.close()
    
    
    
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
    
    print ""
    
    if not skip_listings:
        print "Creating trait listings..."
        traitChi = htmltools.createTraitListingsHTML("results"+os.sep+"html"+os.sep+"traitlists", commonGenes, pfilter_cutoff)
        print "Creating gene listings..."
        htmltools.createGeneListingsHTML("results"+os.sep+"html"+os.sep+"genelists", gwasDB.__geneSet, traitChi)
    
    print "Creating referenced HTML gene and trait list reports..."
    
    # write trait frequency tables
    writeTraitPage("results\\html\\gwas_traits.html", "GWAS Traits for Rapidly Evolving Genes", "Listing of disease traits associated with genes in the rapidly evolving geneset", commonGenes, pfilter_cutoff)
    writeTraitPage("results\\html\\drugbank_traits.html", "GWAS Traits for Drugbank Genes", "Listing of disease traits associated with genes in the drugbank targets database", gwasDB.__geneSet & drugDB.__geneSet, pfilter_cutoff)
    writeTraitPage("results\\html\\overlap_traits.html", "GWAS Traits for overlap of RE and Drugbank", "Listing of disease traits associated with genes in the rapidly evolving geneset and the drugbank targets database", overlap, pfilter_cutoff)
    
    # write gene tables
    
    writeGenePage("results\\html\\gwas_genes.html", "Rapidly Evolving Genes found in GWAS Disease Studies", "Listing of rapidly evolving genes found in GWAS.", commonGenes)
    writeGenePage("results\\html\\drugbank_genes.html", "Rapidly Evolving Genes found in Drugbank", "Listing of rapidly evolving genes found in Drugbank.", commonDrugTargets)
    writeGenePage("results\\html\\drugbank_gwas_genes.html", "Drugbank targets found in GWAS Disease Studies", "Listing of Drugbank targets found in GWAS Disease Studies.", gwas_drugbank_overlap)
    writeGenePage("results\\html\\all_genes.html", "Rapidly Evolving Genes found in GWAS Disease Studies and Drugbank", "Listing of rapidly evolving genes found in both GWAS and Drugbank.", overlap)
    
    
    
    