import os
import sys

gff3_filename=sys.argv[1]
new_gff3_filename=gff3_filename.split(".gff3")[0]+"_modified.gff3"

fhr=open(gff3_filename,"r")
fhw=open(new_gff3_filename,"w")
for line in fhr:
    if "#" in line:
        fhw.write(line)
    else:
        if not ("ID" in line and "Parent" in line):
            fhw.write(line)
            continue
        line=line.strip().split("\t")
        lastcol=line[-1]
        lastcol=lastcol.split(";")
        if lastcol[0].split("ID=")[-1] == lastcol[1].split("Parent=")[-1]:
            lastcol[0]=lastcol[0]+".1"
            lastcol[-1]=lastcol[-1]+".1"
            print("Modifying ",line)
        line[-1]=";".join(lastcol)
        fhw.write("\t".join(line)+"\n")
        
fhw.close()
fhr.close()

cmd="mv "+new_gff3_filename+" "+gff3_filename
#os.system(cmd)