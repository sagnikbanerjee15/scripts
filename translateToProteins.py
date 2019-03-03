"""
Usage translateToProteins.py <name of transcript file> <output filename>
"""

import sys
from Bio.Seq import Seq

inputfilename=sys.argv[1]
outputfilename=sys.argv[2]

fhr=open(inputfilename,"r")
fhw=open(outputfilename,"w")
while True:
    line=fhr.readline()
    if not line :break
    #print(line)
    if len(line.strip())==0:
        continue
    if line[0]==">" and "CDS=" in line:
        id = line.strip().split()[0]
        if "transcript:" in id:
            id=">"+id.split("transcript:")[-1]
        seq=fhr.readline().strip()
        start,end=int(line.strip().split("CDS=")[-1].split("-")[0]),int(line.strip().split("CDS=")[-1].split("-")[-1])
        seq_to_be_translated=seq[start-1:end]
        fhw.write(id+"\n")
        fhw.write(str(Seq(seq_to_be_translated).translate(to_stop=True))+"\n")
fhr.close()
fhw.close()