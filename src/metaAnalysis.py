from geneList import *
import os,sys
import geneVerifier as geneDB
import gwasCatalog as gwasDB



import math
import scipy.stats as stats
import htmltools
from markup import oneliner

__ENABLE_GENE_VERIFICATION = 0
__ENABLE_GENE_UPDATES = 0
__EXCLUDE_PSEUDOGENES = 0
__INCLUDE_MAPPED_GENES = 0
__EXCLUDE_OLFACTORY_PROTEINS = 1
__USE_DRUGBANK_XML = 0


__traitMetaAnalysis = {}

def writeGeneListToFile(filename, genes):
    ffile = open(filename, 'w')
    
    for geneSym in genes:
        ffile.write(geneDB.__original_names[geneSym] + "\n")
    
    ffile.close()

def computeTraitChiSquares(genes, pfilter):
    
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
        
        chisq, oddsratio, kappa = geneUtils.contingentChiSquare(a,b,c,d)
        pvalue = stats.chisqprob(chisq, 1)
        
        traitChi[trait] = (a, chisq, oddsratio, kappa, pvalue, len(traitGenes))
    
    return traitChi
    

def writeGenePage(output_dir, filename, title, desc, total, geneTable):
    genepage = htmltools.createPage(title)

    htmltools.pageDescription(genepage, desc)
    
    newTable=[]
    for row in geneTable:
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
 
def writeTraitPage(filename, title, desc, commonGenes, pfilter_cutoff):
    traitpage = htmltools.createPage(title, scripts={'sorttable.js':'javascript'})
    htmltools.pageDescription(traitpage, desc)
    
    traitChi = computeTraitChiSquares(commonGenes, pfilter_cutoff)
    
    traitTable = []
    
    for trait in traitChi:
        (cnt, chisq, oddsratio, kappa, pvalue,numgenes)=traitChi[trait]
        
        translate = trait.replace(" ","_").replace("/", " or ").replace("\\", " or ")
        
        if len(trait) > 38:
            trait = trait[:35] + "..."
        alink = "<a href=\"traitlists/%s.html\">%s</a>" % (translate, trait)
        
        entry = [alink, cnt, numgenes, "%.2f" % (chisq), "%.7f" % (pvalue), "%.1f" % (oddsratio), "%.4f" % (kappa)]
        traitTable.append(entry)
        
    traitTable = sorted(traitTable, key=lambda item: -item[1])
    
    htmltools.createTable(traitpage, traitTable, ["Disease/Trait", "# RE Genes", "# Trait Genes", "chi-square", "p-value", "odds ratio", "kappa"], 
            "traitlisthead", None, ["traitcol", "recol", "genecol", "chicol",
                "pcol", "oddscol", "kappacol"], "sortable", None)
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
        
        drugs_targeting_RE = set([])
        allDrugs = set([])
        for gene in traitGenes & RE_genes:
            drug_counts_by_gene[gene] = 0
            if gene in drugDB.__drugDict:
                for drug in drugDB.__drugDict[gene]:
                    RE_drugs.append(drug)
                drug_counts_by_gene[gene] += len(drugDB.__drugDict[gene])
                drugs_targeting_RE |= drugDB.__drugDict[gene]

        drugs_targeting_disease = set([])

        for gene in ((traitGenes & drug_genes) - RE_genes):
            drug_counts_by_gene[gene] = 0
            if gene in drugDB.__drugDict:
                for drug in drugDB.__drugDict[gene]:
                    other_drugs.append(drug)
                drug_counts_by_gene[gene] += len(drugDB.__drugDict[gene])
                drugs_targeting_disease |= drugDB.__drugDict[gene]

        drugs_targeting_other_RE = set([])

        for gene in RE_genes:
            drug_counts_by_gene[gene] = 0
            if gene in drugDB.__drugDict:
                drugs_targeting_other_RE |= drugDB.__drugDict[gene]

        drugs_targeting_other_RE -= drugs_targeting_RE

        a = len(drugs_targeting_RE)
        b = len(drugs_targeting_disease)
        c = len(drugs_targeting_other_RE)
        d = len(set(drugDB.__drugs.keys()) - (drugs_targeting_RE |
                drugs_targeting_disease | drugs_targeting_other_RE))
        # print a, b, c, d

        if (a + b) == 0 or (a + c) == 0:
            chi, odds, kappa = 0, 0, 0
            pvalue = 1.0
        else:
            chi, odds, kappa = geneUtils.contingentChiSquare(a,b,c,d)
            pvalue = stats.chisqprob(chi, 1)
        __traitMetaAnalysis[trait]['RE_drugs'] = RE_drugs
        __traitMetaAnalysis[trait]['other_drugs'] = other_drugs
        __traitMetaAnalysis[trait]['drug_counts'] = drug_counts_by_gene
        __traitMetaAnalysis[trait]['drugchi'] = (a,b,c,d,chi,pvalue,odds,kappa)

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
        
        chi_drugs = traitMetadata['drugchi']
        htmltools.createChiTable(traitpage, "Drug contingency for targeting disease vs targeting rapidly evolving proteins:",
                "Targets Disease Genes", "Targets RE Genes", chi_drugs[0], chi_drugs[1],
                chi_drugs[2], chi_drugs[3], chi_drugs[4],
                chi_drugs[5], chi_drugs[6], chi_drugs[7])


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
            else:
                traitpage.li(oneliner.a(drugDB.__drugs[drug]['name'], href=link))
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
            else:
                traitpage.li(oneliner.a(drugDB.__drugs[drug]['name'], href=link))
        traitpage.ul.close()
        traitpage.div.close()


        traitpage.td.close()
        traitpage.tr.close()
        traitpage.table.close()
        
        htmltools.savePage(traitpage, traitListFilename)
    
def createGeneListingsHTML(geneListDir):
    for gene in gwasDB.__geneSet:
        
        genePage = htmltools.createPage("Gene Summary: " +
                geneDB.__original_names[gene], css_file='../genereport.css',
                scripts = {'../sorttable.js':'javascript'})
        
        
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
            oddsratio = __traitMetaAnalysis[trait]['RE_chi'][6]
            kappa     = __traitMetaAnalysis[trait]['RE_chi'][7]
            numgenes  = __traitMetaAnalysis[trait]['geneset_size']
            translate = trait.replace(" ","_").replace("/", " or ").replace("\\", " or ")
            
            if len(trait) > 38:
                trait = trait[:35] + "..."
            traitTable.append(["<a href=\"../traitlists/%s.html\">%s</a>" %
                (translate,trait), cnt, numgenes, "%.2f" % (chisq), "%.7f" %
                (pvalue), "%.1f" % (oddsratio), "%.4f" % (kappa)])
            
        genePage.div("Gene %s, total traits: %d" % (geneDB.__original_names[gene],len(traitTable)), class_="header")
        
        htmltools.createTable(genePage, traitTable, ["Disease/Trait", "#RE Genes", 
            "#Trait Genes", "Chi-square", "p-value", "oddsratio",
            "kappa"], "traitlisthead", None, ["traitcol","recol", "genecol",
                "chicol","pcol","oddscol","kappacol"], "sortable", None)
        
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
                else:
                    genePage.li(oneliner.a(drugDB.__drugs[drug]['name'], href=link))
            genePage.ul.close()
            genePage.div.close()
            
        
        htmltools.savePage(genePage, os.sep.join([geneListDir, gene + ".html"]))

if __name__ == "__main__":
    
    pfilter_cutoff = 0.05
    trait_exclude_file = 0
    gene_frequency_filter=1
    skip_listings=0
    output_dir = os.sep.join(["results", "html"])
    included_studies = [1,2,3,4,5]
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
        elif sys.argv[i] == "--drugbank-xml":
            __USE_DRUGBANK_XML = 1
        elif sys.argv[i] == "--rm-pseudo":
            __EXCLUDE_PSEUDOGENES = 1
        elif sys.argv[i] == "--rm-study":
            for item in sys.argv[i+1].split(","):
                item = item.strip()
                included_studies.remove(int(item))
        elif sys.argv[i] == "--output-dir":
            output_dir = os.sep.join(["results", sys.argv[i+1]])
        elif sys.argv[i].startswith("--"):
            print "Invalid arg:", sys.argv[i]
            sys.exit(1)
    
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
    

    if __USE_DRUGBANK_XML:
        import drugbankDatabase as drugDB
        print "\nLoading drug bank database from XML..."
        drugDB.loadXML(os.sep.join(["data","drugbank","drugbank.xml"]))
        drugDB.mapTargetNames(os.sep.join(["data","drugbank","target_links.csv"]))
    else:
        import drugbankCatalog as drugDB
        print "\nLoading DrugBank drug catalogue..."
        drugDB.initDruglist(os.sep.join(["data","drugbank","drug_links.csv"]))
        print "\nLoading DrugBank drug target catalogue..."
        drugDB.initTargets(os.sep.join(["data","drugbank","target_links.csv"]), os.sep.join(["data","drugbank","all_target_protein.fasta"]),__ENABLE_GENE_VERIFICATION,__ENABLE_GENE_UPDATES)
    
    commonGenes = studyGenes & gwasDB.__geneSet
    cgWithoutGWAS = studyGenes - gwasDB.__geneSet
    gwasWithoutCG = gwasDB.__geneSet - studyGenes
    
    
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
    
    drugsTargetingRE = drugDB.getDrugsTargetingProteinSet(studyGenes)
    drugsTargetingGWAS = drugDB.getDrugsTargetingProteinSet(gwasDB.getDavidBackgroundSet(pfilter_cutoff))
    allDrugs = set(drugDB.__drugs.keys())

    a4 = len(drugsTargetingRE & drugsTargetingGWAS)
    b4 = len(drugsTargetingGWAS - drugsTargetingRE)
    c4 = len(drugsTargetingRE - drugsTargetingGWAS)
    d4 = len(allDrugs - (drugsTargetingGWAS | drugsTargetingRE))
    chisq4, odds4, kappa4 = geneUtils.contingentChiSquare(a4,b4,c4,d4)
    pvalue4 = stats.chisqprob(chisq4, 1)
    
    foldchange_drugs = (float(a4 + c4) / float(a2 + c2)) / (float(b4) / float(a1 + b1))

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
    os.system("copy %s %s" % (os.sep.join(["html", "sorttable.js"]), output_dir ))
    os.system("copy %s %s" % (os.sep.join(["html", "DAVID", str(gene_frequency_filter), "*"]), os.sep.join([output_dir, "DAVID", ""])))
    
    # make index.html
    
    print "Making index.html..."
    
    # Report data set properties
    
    indexpage = htmltools.createPage("Rapidly Evolving Gene Meta-Analysis Report")
    
    indexpage.div("Summary of loaded information:", class_="header")
    
    indexpage.div("Rapidly Evolving genes loaded:                  %d".replace(" ", "&nbsp;") % (len(studyGenes)),class_="console")
    indexpage.div("Drugbank drug target proteins loaded:           %d".replace(" ", "&nbsp;") % (len(drugDB.__geneSet)),class_="console")
    indexpage.div("GWAS Genes indicated in disease studies:        %d".replace(" ", "&nbsp;") % (len(gwasDB.__geneSet)),class_="console")
    
    
    # report gwas, RE overlap chi matrix
    
    htmltools.createChiTable(indexpage, "GWAS Overlap with Rapidly Evolving Geneset:", "GWAS", "RE", a1, b1, c1, d1, chisq1, pvalue1,
            oddsratio1, kappa1)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="gwas_traits.html")
    indexpage.a("Gene Listing", href="gwas_genes.html")
    
    if os.path.exists(os.sep.join([output_dir, "DAVID", "david_gwas_common.xhtml"])):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_gwas_common.xhtml")
        
    indexpage.div.close()
    
    
    # report drugbank, RE overlap chi matrix
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "Drugbank Overlap with Rapidly Evolving Geneset:", "Drugbank", "RE", a2, b2, c2, d2, chisq2, pvalue2, oddsratio2, kappa2)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="overlap_traits.html")
    indexpage.a("Gene Listing", href="drugbank_genes.html")
    
    if os.path.exists(os.sep.join([output_dir, "DAVID", "david_drugbank_common.xhtml"])):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_drugbank_common.xhtml")
    
    indexpage.div.close()
    
    htmltools.createChiTable(indexpage, "Drug contingency for targeting disease vs targeting rapidly evolving proteins:",
            "Targets Disease Genes", "Targets RE Genes", a4, b4, c4, d4,
            chisq4, pvalue4, odds4, kappa4)

    indexpage.div("Fold change drugs per gene (RE) to drugs per gene (GWAS):  %.2f" % (foldchange_drugs), class_ = "description")
    # report drugbank, GWAS overlap chi matrix
    
    indexpage.div.open(class_="reportsquare")
    htmltools.createChiTable(indexpage, "GWAS Overlap with Drugbank:",
            "Drugbank", "GWAS", a3, b3, c3, d3, chisq3, pvalue3, 
            oddsratio3, kappa3)
    
    indexpage.div.open(class_="links")
    indexpage.a("Trait Report", href="drugbank_traits.html")
    indexpage.a("Gene Listing", href="drugbank_gwas_genes.html")
    
    if os.path.exists(os.sep.join([output_dir, "DAVID", "david_drugbank_gwas_common.xhtml"])):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_drugbank_gwas_common.xhtml")
    
    indexpage.div.close()
    
    
    
    
    overlap = commonGenes & studyGenes & drugDB.__geneSet
    
    indexpage.div.open(class_="links")
    indexpage.p.open()
    indexpage.add("Total overlap of drugbank, GWAS, and Rapidly Evolving geneset: ")
    indexpage.a(str(len(overlap)), href="all_genes.html")
    
    if os.path.exists(os.sep.join([output_dir, "DAVID", "david_all.xhtml"])):
        indexpage.br()
        indexpage.a("DAVID Results", href="DAVID/david_all.xhtml")
    
    indexpage.p.close()
    indexpage.div.close()
    
    htmltools.endPage(indexpage)
    htmltools.savePage(indexpage, os.sep.join([output_dir, "index.html"]))
    
    print ""
    
    if not skip_listings:
        print "Creating trait listings..."
        computeTraitGeneLists(studyGenes, drugDB.__geneSet, pfilter_cutoff)
        computeTraitDrugLists(studyGenes, drugDB.__geneSet, pfilter_cutoff)
        createTraitListingsHTML(os.sep.join([output_dir, "traitlists"]))
        print "Creating gene listings..."
        createGeneListingsHTML(os.sep.join([output_dir,"genelists"]))
    
    print "Creating referenced HTML gene and trait list reports..."
    
    # write trait frequency tables
    writeTraitPage(os.sep.join([output_dir,"gwas_traits.html"]), "GWAS Traits for Rapidly Evolving Genes", "Listing of disease traits associated with genes in the rapidly evolving geneset", studyGenes, pfilter_cutoff)
    writeTraitPage(os.sep.join([output_dir,"drugbank_traits.html"]), "GWAS Traits for Drugbank Genes", "Listing of disease traits associated with genes in the drugbank targets database", gwasDB.__geneSet & drugDB.__geneSet, pfilter_cutoff)
    writeTraitPage(os.sep.join([output_dir,"overlap_traits.html"]), "GWAS Traits for overlap of RE and Drugbank", "Listing of disease traits associated with genes in the rapidly evolving geneset and the drugbank targets database", overlap, pfilter_cutoff)
    
    # write gene tables
    
    total, geneTable = getGWASFrequencyTable(commonGenes)
    writeGenePage(output_dir, os.sep.join([output_dir,"gwas_genes.html"]), "Rapidly Evolving Genes found in GWAS Disease Studies", "Listing of rapidly evolving genes found in GWAS.", total, geneTable)
    writeGeneListToFile(os.sep.join(["results","log","gwas_genes.txt"]), [row[0] for row in geneTable])

    total, geneTable = getGWASFrequencyTable(commonDrugTargets)
    writeGenePage(output_dir, os.sep.join([output_dir,"drugbank_genes.html"]), "Rapidly Evolving Genes found in Drugbank", "Listing of rapidly evolving genes found in Drugbank.", total, geneTable)
    writeGeneListToFile(os.sep.join(["results","log","drugbank_genes.txt"]), [row[0] for row in geneTable])
    
    total, geneTable = getGWASFrequencyTable(gwas_drugbank_overlap)
    writeGenePage(output_dir, os.sep.join([output_dir,"drugbank_gwas_genes.html"]), "Drugbank targets found in GWAS Disease Studies", "Listing of Drugbank targets found in GWAS Disease Studies.", total, geneTable)
    writeGeneListToFile(os.sep.join(["results","log","drugbank_gwas_genes.txt"]), [row[0] for row in geneTable])
    
    total, geneTable = getGWASFrequencyTable(overlap)
    writeGenePage(output_dir, os.sep.join([output_dir,"all_genes.html"]), "Rapidly Evolving Genes found in GWAS Disease Studies and Drugbank", "Listing of rapidly evolving genes found in both GWAS and Drugbank.", total, geneTable)
    writeGeneListToFile(os.sep.join(["results","log","all_overlap_genes.txt"]), [row[0] for row in geneTable])
    
    
    
    
