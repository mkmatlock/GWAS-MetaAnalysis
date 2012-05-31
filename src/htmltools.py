import markup
import geneVerifier as geneDB
import gwasCatalog as gwasDB
import drugbankCatalog as drugDB
import os
import geneUtils
import scipy.stats as stats
import geneVerifier as geneDB




    
def savePage( page, filename ):
    ofile = open(filename,'w')
    ofile.write(str(page))
    ofile.close()
    
def createPage( pagetitle,css_file='genereport.css',scripts={}):
    page = markup.page()
    page.init( title=pagetitle,
               css = (css_file),
               script=scripts)

    page.div.open(id="body")
    
    return page

def endPage(page):
    page.div.close()
               
def pageDescription( page, desc ):
    page.div("Description:", id="dheader")
    page.div(desc, id="description")
               
def createChiTable(page, name, cat1, cat2, a, b, c, d, chisq, pvalue, oddsratio, kappa, fisherp):
    page.div(name,class_="header")
    page.div.open(class_="chireport")
    page.p("Test Matrix:")
    createTable(page, [["in "+cat1, a, b],["not in "+cat1, c, d]], ["","in "+cat2,"not in "+cat2], "chiheader", ["chirow1","chirow2"], ["chicol1","chicol2","chicol3"], None, "chimatrix")
    createTable(page, [["Chi statistic", "%.2f" % (chisq)],["P-value", "%.7f" %
        (pvalue)], ["Odds-Ratio", "%.1f" % (oddsratio)], ["Kappa Statistic",
            "%.4f" % (kappa)], ["Fisher P-value", "%.7f" % (fisherp)]], col_classes=["name", "value"], table_cls = "chireportstats")
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
