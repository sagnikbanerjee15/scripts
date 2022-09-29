"""
Program will find structures which are on incorrect 
strands
"""
import os
import sys

gff3_filename=sys.argv[1]
new_gff3_filename=gff3_filename.split(".gff3")[0]+"_modified.gff3"

fhr=open(gff3_filename,"r")
parent_to_strand={}
all_parents=[]
for line in fhr:
    if "#" in line:continue
    line=line.strip().split("\t")
    #print(line)
    strand=line[6]
    lastcol=line[-1]
    if "Parent" in lastcol:
        parent_id=lastcol.split("Parent=")[-1].split(";")[0]
        if parent_id+"_"+strand not in parent_to_strand:
            parent_to_strand[parent_id+"_"+strand]=0
        parent_to_strand[parent_id+"_"+strand]+=1
        all_parents.append(parent_id)
fhr.close()

parents_to_be_changed={}
for parent in list(set(all_parents)):
    if parent+"_+"  in parent_to_strand and parent+"_-" in parent_to_strand:
        print(parent,parent_to_strand[parent+"_+"],parent_to_strand[parent+"_-"])
        parents_to_be_changed[parent]="+" if parent_to_strand[parent+"_+"] >= parent_to_strand[parent+"_-"] else "-"

fhr=open(gff3_filename,"r")
fhw=open(new_gff3_filename,"w")
for line in fhr:
    if "#" in line:
        fhw.write(line)
        continue
    if line.strip().split("\t")[-1].split("Parent=")[-1].split(";")[0] not in parents_to_be_changed:
        fhw.write(line)
    else:
        line=line.strip().split("\t")
        line[6]=parents_to_be_changed[line[-1].split("Parent=")[-1].split(";")[0]]
        fhw.write("\t".join(line)+"\n")
fhw.close()
fhr.close()
