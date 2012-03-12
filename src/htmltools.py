import markup
import geneVerifier as geneDB
import gwasCatalog as gwasDB
import drugbankCatalog as drugDB
import os
import geneUtils
import scipy.stats as stats
import geneVerifier as geneDB


def createGeneListTable(page, genes, verify=None):
    output_list = []
    
    total = 0
    for geneSym in genes:
        l = len(gwasDB.__traitDict[geneSym])
        total += l
        
        if verify==None or geneSym in verify:
            output_list.append(["<a href=\"../genelists/%s.html\" > %s </a>" % (geneSym,geneDB.__original_names[geneSym]), l])
        else:
            output_list.append([geneDB.__original_names[geneSym], l])
        
    output_list = sorted(output_list, key=lambda item: -item[1])
    createTable(page, output_list, ["Gene","#Associated Traits"], "genelisthead", None, ["genecol","traitcol"],"genetable",None)

    
def savePage( page, filename ):
    ofile = open(filename,'w')
    ofile.write(str(page))
    ofile.close()
    
def createPage( pagetitle,css_file='../genereport.css'):
    page = markup.page()
    page.init( title=pagetitle,
               css = (css_file))

    page.div.open(id="body")
    
    return page

def endPage(page):
    page.div.close()
               
def pageDescription( page, desc ):
    page.div("Description:", id="dheader")
    page.div(desc, id="description")
               
def createChiTable(page, name, cat1, cat2, a, b, c, d, chisq, pvalue):
    page.div(name,class_="header")
    page.div.open(class_="chireport")
    page.p("Test Matrix:")
    createTable(page, [["in "+cat1, a, b],["not in "+cat1, c, d]], ["","in "+cat2,"not in "+cat2], "chiheader", ["chirow1","chirow2"], ["chicol1","chicol2","chicol3"], None, "chimatrix")
    createTable(page, [["Chi statistic", "%.2f" % (chisq)],["P-value", "%.7f" % (pvalue)]], col_classes=["name", "value"], table_cls = "chireportstats")
    page.div.close()
    
    
    
def createTable( page, data, header = None, header_class = None, row_classes=None, col_classes=None, table_cls = None, table_id = None ):
    
    if table_cls != None and table_id != None:
        page.table.open( class_=table_cls, id=table_id)
    elif table_cls != None:
        page.table.open( class_=table_cls )
    elif table_id != None:
        page.table.open( id=table_id )
    else:
        page.table.open( )
        
    if header != None:
        page.tr( class_=header_class )
        for item in header:
            page.th( item )
        page.tr.close( )
    
    i=0
    for row in data:
        if row_classes != None:
            page.tr( class_=row_classes[i] )
            i+=1
        else:
            page.tr.open( )
        
        j=0
        for col in row:
            if col_classes != None:
                page.td( col, class_=col_classes[j] )
                j+=1
            else:
                page.td( col )
        
        page.tr.close( )
    
    page.table.close( )