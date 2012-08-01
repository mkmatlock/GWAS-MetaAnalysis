__DEBUG=0

def parseFASTA(fasta_file):
    fasta_file = open(fasta_file,'r')
    
    seq = []
    cur = 0
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):
            type = line[line.find(">")+1 : line.find("|")]
            description = line[line.find("|")+1:]
            cur = [type,description,""]
            seq.append(cur)
        elif cur!=0:
            cur[2] += line
    
    fasta_file.close()
    return seq