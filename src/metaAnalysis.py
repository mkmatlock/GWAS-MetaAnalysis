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
    
    
    # chi-square contingency table for gwas and RE
    
    a1 = len(commonGenes)
    b1 = len(gwasWithoutCG)
    c1 = len(cgWithoutGWAS)
    d1 = len( geneDB.__approved_symbols - (studyGenes | gwasDB.__geneSet) )
    chisq1 = geneUtils.contingentChiSquare(a,b,c,d)
    pvalue1 = stats.chisqprob(chisq, 1)
    
    
    a2 = len( studyGenes & drugDB.__geneSet )
    b2 = len( drugDB.__geneSet - studyGenes )
    c2 = len( studyGenes - drugDB.__geneSet )
    d2 = len( geneDB.__approved_symbols - ( studyGenes | drugDB.__geneSet ) )
    chisq2  = geneUtils.contingentChiSquare(a,b,c,d)
    pvalue2 = stats.chisqprob(chisq, 1)
    
    
    
    commonDrugTargets = drugDB.__geneSet & studyGenes
    gwas_drugbank_overlap=gwasDB.__geneSet & drugDB.__geneSet
    
    a3 = len( gwas_drugbank_overlap )
    b3 = len( drugDB.__geneSet - gwasDB.__geneSet )
    c3 = len( gwasDB.__geneSet - drugDB.__geneSet )
    d3 = len( geneDB.__approved_symbols - ( drugDB.__geneSet | gwasDB.__geneSet ) )
    chisq3  = geneUtils.contingentChiSquare(a,b,c,d)
    pvalue3 = stats.chisqprob(chisq, 1)
    
    
    # start preparing HTML reports
    
    if not os.path.exists("results" + os.sep + "html"):
        os.mkdir("results" + os.sep + "html")
        os.mkdir("results" + os.sep + "html" + os.sep + "genelists")
        os.mkdir("results" + os.sep + "html" + os.sep + "traitlists")
        os.mkdir("results" + os.sep + "html" + os.sep + "DAVID")
        os.system("cp %s %s" % ("html" + os.sep + "genereport.css", "results" + os.sep + "html")
        os.system("cp %s %s" % ("html" + os.sep + "DAVID_" + gene_frequency_filter + os.sep + "*", "results" + os.sep + "html" + os.sep + "DAVID")
    
    # make index.html
    
    # Report data set properties
    
    indexpage = htmltools.createPage("Rapidly Evolving Gene Meta-Analysis Report")
    
    indexpage.div("Summary of loaded information:", class_="header")
    
    indexpage.div("Rapidly Evolving genes loaded:                  %d".replace(" ", "&nbsp;") % (len(studyGenes)),class_="console")
    indexpage.div("Drugbank drug target proteins loaded:           %d".replace(" ", "&nbsp;") % (len(drugDB.__geneSet)),class_="console")
    indexpage.div("GWAS Genes indicated in disease studies:        %d".replace(" ", "&nbsp;") % (len(gwasDB.__geneSet)),class_="console")
    
    
    # report gwas, RE overlap chi matrix
    
    htmltools.createChiTable(indexpage, "GWAS Overlap with Rapidly Evolving Geneset:", "GWAS", "RE", a1, b1, c1, d1, chisq1, pvalue1)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="gwas_traits.html")
    indexpage.a("Gene Listing", href="gwas_genes.html")
    
    if os.path.exists("results" + os.sep + "html" + os.sep + "DAVID" + os.sep + "david_gwas_common.xhtml"):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_gwas_common.xhtml")
        
    indexpage.div.close()
    indexpage.div.close()
    
    
    # report drugbank, RE overlap chi matrix
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "Drugbank Overlap with Rapidly Evolving Geneset:", "Drugbank", "RE", a2, b2, c2, d2, chisq2, pvalue2)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="overlap_traits.html")
    indexpage.a("Gene Listing", href="drugbank_genes.html")
    
    if os.path.exists("results" + os.sep + "html" + os.sep + "DAVID" + os.sep + "david_drugbank_common.xhtml"):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_drugbank_common.xhtml")
    
    indexpage.div.close()
    indexpage.div.close()
    
    
    # report drugbank, GWAS overlap chi matrix
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "GWAS Overlap with Drugbank:", "Drugbank", "GWAS", a3, b3, c3, d3, chisq3, pvalue3)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="drugbank_traits.html")
    indexpage.a("Gene Listing", href="drugbank_gwas_genes.html")
    
    if os.path.exists("results" + os.sep + "html" + os.sep + "DAVID" + os.sep + "david_drugbank_gwas_common.xhtml"):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_drugbank_gwas_common.xhtml")
    
    indexpage.div.close()
    indexpage.div.close()
    
    
    overlap = commonGenes & studyGenes & drugDB.__geneSet
    
    indexpage.div.open(class_="links")
    indexpage.p.open()
    indexpage.add("Total overlap of drugbank, GWAS, and Rapidly Evolving geneset: ")
    indexpage.a(str(len(overlap)), href="all_genes.html")
    
    if os.path.exists("results" + os.sep + "html" + os.sep + "DAVID" + os.sep + "david_all.xhtml"):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_all.xhtml")
    
    indexpage.p.close()
    indexpage.div.close()
    
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
    
    
    
    