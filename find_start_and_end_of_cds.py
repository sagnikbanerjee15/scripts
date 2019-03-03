import argparse
import sys
import os
import re

def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python find_start_and_end_of_cds.py --help for more details
    """
    parser = argparse.ArgumentParser(prog="find_start_and_end_of_cds.py",description="Find the start and the end coordinates of transcripts which code a protein. The output will be in the form of a csv file")
    parser.add_argument("--gtf","-gtf",help="Enter the gtf file",required=True)
    parser.add_argument("--cdna","-cdna",help="Enter the filename of the cDNA file",required=True)
    parser.add_argument("--cds","-cds",help="Enter the filename of the CDS file",required=True)
    parser.add_argument("--functions","-f",help="Enter the filename with transcript to function mappings",required=True)
    return parser.parse_args()

def extractTranscriptToFunctionMappings(filename):
    """
    Returns a dictionary with transcript to 
    function mappings
    """
    transcript_to_function={}
    fhr=open(filename,"r")
    for line in fhr:
        transcript,function=line.strip().split(",")
        transcript_to_function[transcript]=function
    fhr.close()
    return transcript_to_function

def readFastaFile(filename):
    """
    Reads in the fasta file and returns a dictionary
    """  
    data={}
    fhr=open(filename,"r")
    for line in fhr:
        if ">" in line:
            data[line[1:].split()[0]]=fhr.readline().strip()
    fhr.close()
    return data

def locateStartAndStopOfEachTranscript(cdna,cds):
    """
    Returns 1-based start and end coordinates of 
    the protein coding transcript
    """
    return [(m.start(0)+1, m.end(0)) for m in re.finditer(cds, cdna)]

def extractTranscriptToGeneMapsFromGTFFile(filename):
    """
    Generate mappings from transcript to gene
    """
    transcript_to_gene_map={}
    fhr=open(filename,"r")
    for line in fhr:
        if line[0]=="#":continue
        chromosome=line.strip().split("\t")[0]
        info=line.strip().split("\t")[-1].split(";")
        #print(info)
        transcript,gene=info[0],info[1]
        #print(gene,transcript)
        transcript=transcript.split()[-1].strip("\"")
        gene=gene.split()[-1].strip("\"")
        transcript_to_gene_map[transcript]=[gene,chromosome]
        #print(chromosome,gene,transcript)
    fhr.close()
    return transcript_to_gene_map
    
def findStartAndStopCoordinates(options):
    """
    
    """
    transcript_to_gene=extractTranscriptToGeneMapsFromGTFFile(options.gtf)
    cdna=readFastaFile(options.cdna)
    cds=readFastaFile(options.cds)
    transcript_to_function=extractTranscriptToFunctionMappings(options.functions)
    #print(list(transcript_to_gene.keys()))
    for transcript in cdna:
        #if transcript!="AGP50780.1":continue
        if transcript in cds:
            #print(transcript_to_gene[transcript])
            start,end=list(locateStartAndStopOfEachTranscript(cdna[transcript], cds[transcript])[0])
            print(",".join(map(str,[transcript_to_gene[transcript][1],transcript,transcript_to_gene[transcript][0],start,end,len(cdna[transcript]),transcript_to_function[transcript] if transcript in transcript_to_function else ""])))
            
def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    findStartAndStopCoordinates(options)

if __name__ == "__main__":
    main()