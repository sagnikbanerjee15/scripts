import os
import argparse
import sys
import math
import gc
import numpy as np

def parseCommandLineArguments():
    """
    """
    parser = argparse.ArgumentParser(prog="convert_fastq_to_fasta.py",description="Convert fastq to fasta file")
    parser.add_argument("--input","-i",help="Enter the full path of the fastq file",required=True)
    parser.add_argument("--prefix","-p",help="Enter the prefix of the output file",required=True)
    parser.add_argument("--unique","-u",nargs='*',default=False)
    return parser.parse_args()

def readFastqFile(filename):
    """
    Reads in the fastq file and returns a dictionary object
    """
    data={}
    fhr=open(filename)
    cntr=0
    while(True):
        line=fhr.readline().strip()
        if not line: break
        #if cntr%4==0:
        id = line.strip()[1:]
        seq = fhr.readline().strip()
        useless=fhr.readline().strip()
        quality = fhr.readline().strip()
        data[id]={"seq":seq,"useless":useless,"quality":quality}
    #print(list(data.keys())[:5])
    return data

def convertFastqToFasta(data,filename):
    """
    Takes a dictionary as input and write the reads to the file
    """
    fhw=open(filename+".fasta","w")
    for id in data:
        fhw.write(">"+data[id]["seq"]+"\n")
    fhw.close()

def findUniqueReads(data,prefix):
    """
    Looks into the data and flag the reads which have duplicates.
    """
    sequences={}
    sequences_with_quality={}
    for id in data:
        if data[id]["seq"] not in sequences:
            sequences[data[id]["seq"]]=id
        if data[id]["seq"]+"_"+data[id]["quality"] not in sequences_with_quality:
            sequences_with_quality[data[id]["seq"]+"_"+data[id]["quality"]]=id
    
    fhw=open(prefix+"_unique.fasta","w")
    for seq in sequences:
        fhw.write(">"+sequences[seq]+"\n"+seq+"\n")
    fhw.close()
    
    fhw=open(prefix+"_unique.fastq","w")
    for seq in sequences_with_quality:
        fhw.write("@"+sequences_with_quality[seq]+"\n"+seq.split("_")[0]+"\n+\n"+seq.split("_")[1]+"\n")
    fhw.close()

def computeAverageQuality(quality):
    """
    Returns the average quality
    """
    total=0
    for ele in quality:
        total+=ord(ele)
    return total/len(quality)

def retainBest(inputfilename,prefix):
    """
    Retains the best quality alignment among all those having the same sequence
    """
    fhr=open(inputfilename)
    sequences_with_quality={}
    seq_processed=0
    while(True):
        line=fhr.readline().strip()
        if not line: break
        #if cntr%4==0:
        id = line.strip()[1:]
        seq = fhr.readline().strip()
        useless=fhr.readline().strip()
        quality = fhr.readline().strip()
        seq_processed+=1
        #data[id]={"seq":seq,"useless":useless,"quality":quality}
        if seq not in sequences_with_quality:
            sequences_with_quality[seq]={"id":id,"avg_qual":computeAverageQuality(quality),"quality":quality}
        else:
            if computeAverageQuality(quality) > sequences_with_quality[seq]["avg_qual"]:
                sequences_with_quality[seq]={"id":id,"avg_qual":computeAverageQuality(quality),"quality":quality}
        if len(sequences_with_quality)%10000==0:
            print(len(sequences_with_quality),seq_processed,len(sequences_with_quality)/seq_processed)
            sys.stdout.flush()
        """if len(sequences_with_quality)==10000000:
            break"""
    try:
        #fhw=open(prefix+"_best.fastq","w")
        with open(prefix+"_best.fastq","w") as fhw:
            for num,seq in enumerate(list(sequences_with_quality.keys())):
                fhw.write("@"+sequences_with_quality[seq]["id"]+"\n"+seq+"\n+\n"+sequences_with_quality[seq]["quality"]+"\n")
                #del sequences_with_quality[seq]
                if num%1000==0:
                    fhw.flush()
                    #gc.collect()
    finally:
        fhw.close()
        
def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    #data=readFastqFile(options.input)
    data=[]
    if options.unique==False:
        convertFastqToFasta(data, options.prefix)
    else:
        #findUniqueReads(data,options.prefix)
        retainBest(options.input,options.prefix)
    

if __name__ == "__main__":
    main()