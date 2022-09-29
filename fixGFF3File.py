"""
Program fixes gff3 files by replacing repeated IDs with unique ones.
Also removes words like gene: and transcript: from the file
"""
import sys

inputfilename=sys.argv[1]
outputfilename=sys.argv[2]

fhr=open(inputfilename,"r")
fhw=open(outputfilename,"w")
for line in fhr:
    if "#" in line:
        fhw.write(line)
    if "ID" in line:
        line.strip.split("\t")[-1]

fhw.close()
fhr.close()