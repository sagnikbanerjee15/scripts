import sys
import os

gff3_filename=sys.argv[1]
new_gff3_filename=gff3_filename+".temp"
cdna_filename=sys.argv[2]
new_cdna_filename=cdna_filename+".temp"
cds_filename=sys.argv[3]
new_cds_filename=cds_filename+".temp"
prot_filename=sys.argv[4]
new_prot_filename=prot_filename+".temp"

fhr=open(cdna_filename,"r")
fhw=open(new_cdna_filename,"w")
for line in fhr:
    if "Pt" in line:
        newline=line.split()[0]+".1 "+" ".join(line.split()[1:])
        fhw.write(newline.strip()+"\n")
    else:
        fhw.write(line.strip()+"\n")
fhr.close()
fhw.close()

fhr=open(cds_filename,"r")
fhw=open(new_cds_filename,"w")
for line in fhr:
    if "Pt" in line:
        newline=line.split()[0]+".1 "+" ".join(line.split()[1:])
        fhw.write(newline.strip()+"\n")
    else:
        fhw.write(line.strip()+"\n")
fhr.close()
fhw.close()

fhr=open(prot_filename,"r")
fhw=open(new_prot_filename,"w")
for line in fhr:
    if "Pt" in line:
        newline=line.split()[0]+".1 "+" ".join(line.split()[1:])
        fhw.write(newline.strip()+"\n")
    else:
        fhw.write(line.strip()+"\n")
fhr.close()
fhw.close()

fhr=open(gff3_filename,"r")
fhw=open(new_gff3_filename,"w")
for line in fhr:
    if line.strip().split()[0]!="Pt":
        fhw.write(line)
    else:
        if "RNA" in line.strip().split("\t")[2] and "gene" not in line.strip().split("\t")[2]:
            line=line.strip().split("\t")
            last_col=line[-1]
            last_col=last_col.split(";")
            last_col[0]=last_col[0]+".1"
            last_col[-1]=last_col[-1]+".1"
            last_col=";".join(last_col)
            line[-1]=last_col
            fhw.write("\t".join(line)+"\n")
        elif line.strip().split("\t")[2]=="CDS":
            line=line.strip().split("\t")
            last_col=line[-1]
            last_col=last_col.split(";")
            last_col[1]=last_col[1]+".1"
            #last_col[-1]=last_col[-1]+".1"
            last_col=";".join(last_col)
            line[-1]=last_col
            fhw.write("\t".join(line)+"\n")
        elif line.strip().split("\t")[2]=="exon":
            line=line.strip().split("\t")
            last_col=line[-1]
            last_col=last_col.split(";")
            last_col[0]=last_col[0]+".1"
            #last_col[-1]=last_col[-1]+".1"
            last_col=";".join(last_col)
            line[-1]=last_col
            fhw.write("\t".join(line)+"\n")
        else:
            fhw.write(line)
fhw.close()
fhr.close()

cmd="mv "+new_gff3_filename+" "+gff3_filename
os.system(cmd)

cmd="mv "+new_cds_filename+" "+cds_filename
os.system(cmd)

cmd="mv "+new_cdna_filename+" "+cdna_filename
os.system(cmd)

cmd="mv "+new_prot_filename+" "+prot_filename
os.system(cmd)