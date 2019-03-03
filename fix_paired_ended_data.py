"""
Fixes paired ended libraries by associating reads on correct lines
"""

import argparse
import os
import sys

def readBothFastqFiles(filename1,filename2):
    """
    Reads in both the files
    """
    whole_data={}
    fhr1=open(filename1,"r")
    while True:
        id=fhr1.readline().strip()
        if not id:
            break
        seq=fhr1.readline().strip()
        useless=fhr1.readline().strip()
        qual=fhr1.readline().strip()
        whole_data[id]={}
        whole_data[id]["firstpair"]=[seq,useless,qual]
    fhr1.close()
    
    fhr2=open(filename2,"r")
    while True:
        id=fhr2.readline().strip()
        if not id:
            break
        seq=fhr2.readline().strip()
        useless=fhr2.readline().strip()
        qual=fhr2.readline().strip()
        if id not in whole_data:
            whole_data[id]={}
        whole_data[id]["secondpair"]=[seq,useless,qual]
        #whole_data[id]={"secondpair":[seq,useless,qual]}
    fhr2.close()
    return whole_data

def writeBothFastqFiles(whole_data,prefix):
    """
    """
    discarded=[]
    fhw1=open(prefix+"_1.fastq","w")
    fhw2=open(prefix+"_2.fastq","w")
    for id in whole_data:
        #print(whole_data[id])
        if "firstpair" not in whole_data[id] or "secondpair" not in whole_data[id]:
            discarded.append(id)
            continue
        fhw1.write(id+"\n"+whole_data[id]["firstpair"][0]+"\n"+whole_data[id]["firstpair"][1]+"\n"+whole_data[id]["firstpair"][2]+"\n")
        fhw2.write(id+"\n"+whole_data[id]["secondpair"][0]+"\n"+whole_data[id]["secondpair"][1]+"\n"+whole_data[id]["secondpair"][2]+"\n")
    fhw1.close()
    fhw2.close()
    return discarded

def parseCommandLineArguments():
    """
    """
    parser = argparse.ArgumentParser(prog="fix_paired_ended_data.py",description="Fixes paired ended data by correcting the order of the reads in the which they appear. The ids of both the reads must be same for this to work.")
    parser.add_argument("--read1","-r1",help="Please enter the filename of the first read pairs",required=True)
    parser.add_argument("--read2","-r2",help="Please enter the filename of the second read pairs",required=True)
    parser.add_argument("--output_prefix","-out",help="Please enter the output prefix",required=True)
    return parser.parse_args()

def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    whole_data=readBothFastqFiles(options.read1,options.read2)
    discarded=writeBothFastqFiles(whole_data,options.output_prefix)
    print(len(discarded),"reads were discarded")

if __name__ == "__main__":
    main()