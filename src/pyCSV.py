import gc
import os
import re
import codecs

class pyCSV:
    def __init__(self):
        self.rows = 0
        self.cols = 0
        
        self.struct = {}
        self.rowNames = {}
        self.colNames = {}
        
        self.samplerate = 1000
        
        self.name = 0
        
    def set(self, row, col, val):
        if not self.struct.has_key(row):
            self.struct[row] = {}
            
            if self.rows < row:
                self.rows = row
            
        self.struct[row][col] = val
        
        if self.cols < col:
            self.cols = col
    
    def get(self, row, col):
        if not self.struct.has_key(row):
            return None
        if not self.struct[row].has_key(col):
            return None
        
        return self.struct[row][col]
    
    def setRowName(self, row, name):
        self.rowNames[row] = name
    
    def setColName(self, col, name):
        self.colNames[col] = name  
    
    def getRowName(self, row):
        if row not in self.rowNames:
            return None
        return self.rowNames[row]
    
    def getColName(self, col):
        if col not in self.colNames:
            return None
        return self.colNames[col]
    
    def getColumn(self, col):
        column_data = []
        for i in xrange(srow, self.rows+1):
            column_data.append(self.get(i, col))
            
        return column_data
        
    def replaceColumn(self, col, data):
        for i in xrange(0, self.rows+1):
            self.set(i, col, data[i])
    
    def applyFunction(self, col, func):
        for row in self.struct:
            items = self.struct[row]
            if items.has_key(col):
                items[col] = func(items[col])
    
    def downSample(self, ratio):
        newstruct = {}
        newrownames = {}
        for i in xrange(0, (self.rows+1) / ratio):
            row = i * ratio
            
            if self.rowNames.has_key(row):
                newrownames[i] = self.rowNames[row]
            
            newstruct[i] = {}
            
            for j in xrange(0, self.cols+1):
                newstruct[i][j] = self.get(row,j)
                
            del self.struct[row]
        
        del self.struct
        
        self.struct = newstruct
        self.rowNames = newrownames
        self.samplerate /= ratio
        self.rows = self.rows / ratio
    
    def getAsCSVString(self, row, col):
        item = self.get(row,col)
        
        str_r = ""
        if item != None:
            str_r = str(item)
            if(str_r.find(",") > -1):
                str_r = "\"" + str_r + "\""
        return str_r
    
    def writeItem(self, file, str_r):
        file.write( str( str_r ) )
        self.wiPtr = self.writeItemComma
        
    def writeItemComma(self,file,str_r):
        file.write( "," + str(str_r) )
    
    def save(self, filename):
        file = open(filename, 'w')
        
        row = 0
        for i in xrange(0, self.rows+1):
            self.wiPtr = self.writeItem
            for j in xrange(0, self.cols+1):
                self.wiPtr(file, self.getAsCSVString(i,j))
                
            file.write("\n")
        
        file.close()
    
    def transferDownsample(self, input, output, ratio):
        ifile = open(input, 'r')
        ofile = open(output, 'w')
        
        r=0
        report = 100000
        for line in ifile:
            
            if r % ratio == 0:
                ofile.write(line)
            
            if r % report == 0:
                print "."
            
            r+=1
        ifile.close()
        ofile.close()
    
    def load(self, filename, delimiter = ",", codec = 0):
    
        file = 0
        if codec != 0:
            file = codecs.open(filename, 'r', codec)
        else:
            file = open(filename, 'r')
        
        self.name = filename[filename.rfind(os.sep)+1:]
        
        row = 0
        for line in file:
            line = line.strip()
            if line.endswith("\n"):
                line = line[:-1]
            items = line.split(delimiter)
            
            rappend = 0
            i = 0
            
            for item in items:
                stripped = item.strip()
                
                if rappend==0:
                    if stripped.startswith("\"") and stripped.endswith("\""):
                        stripped = stripped[1:-1]
                        self.set(row, i, stripped)
                        i+=1
                    elif stripped.startswith("\""):
                        rappend = stripped
                    else:
                        self.set(row, i, stripped)
                        i+=1
                else:
                    rappend += stripped
                    if stripped.endswith("\""):
                        self.set(row, i, rappend[:-1])
                        i+=1
                        rappend = 0
                
            del line
            del items
            
            row += 1
        file.close()
        gc.collect()