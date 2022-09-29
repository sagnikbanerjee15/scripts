"""
Splits each chromosome into smaller chromosomes such the length of the 
sub-chromosomes do not exceed 2^29-1. This is done mainly to view the 
gene annotations on IGV
Requires the following tools
    - samtools
    - genometools 
    - cufflinks
"""

import argparse
import os
import sys
import math

def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python compute_y2h_scores.py --help for more details
    """
    parser = argparse.ArgumentParser(prog="split_large_chromosomes.py",description="Split each chromosome into smaller chunks")
    parser.add_argument("--genome","-g",help="Please enter the name of the genome fasta file. Please make sure the genome file is a singleline fasta.",required=True)
    parser.add_argument("--gene_annotation","-ga",help="Please enter the name of the gene annotation file in gff3 format",required=True)
    parser.add_argument("--out_prefix","-p",help="Please enter the prefix for output",required=True)
    return parser.parse_args()

def findChromosomes(options):
    """
    Return a dictionary with the chromosome information 
    which has to be split
    """
    cmd="samtools faidx "+options.genome
    os.system(cmd)
    
    chromosome_to_length={}
    fai_filename=options.genome+".fai"
    fhr=open(fai_filename,"r")
    for line in fhr:
        chromosome,length=line.strip().split()[0],int(line.strip().split()[1])
        if length >= math.pow(2,28)-1:
            chromosome_to_length[chromosome]=length
    fhr.close()
    return chromosome_to_length

def findRegionsOfNsEachChromosome(sequence):
    """
    """
    regions_with_N=[]
    start=-1
    #print(sequence[:100])
    i=int(math.pow(2,28))
    while i<len(sequence):
        if sequence[i]=="N":
            if start==-1:
                start=i
        if start!=-1 and sequence[i]!="N":
            regions_with_N.append([start,i-1])
            start=-1
        i+=1
    return regions_with_N

def findRegionsOfNs(options,chromosome_to_length):
    """
    Returns a dictionary of chromosomes to regions which are Ns
    """
    chromosome_to_regions_with_N={}
    fhr=open(options.genome,"r")
    for line in fhr:
        if ">" in line:
            chromosome=line.strip().split()[0][1:]
            if chromosome not in chromosome_to_length:
                fhr.readline()
                continue
        sequence=line.strip()
        regions_with_N=findRegionsOfNsEachChromosome(sequence)
        chromosome_to_regions_with_N[chromosome]=regions_with_N
    fhr.close()
    return chromosome_to_regions_with_N

def modifyBedFile(filename):
    """
    Removes lines with chromosomes
    """
    fhw=open(filename.split(".bed")[0]+".temp","w")
    fhr=open(filename,"r")
    for line in fhr:
        if "chromosome" in line:
            if line.strip().split()[7]=="chromosome":continue
        fhw.write(line)
    fhr.close()
    fhw.close()
    cmd="mv "+filename.split(".bed")[0]+".temp "+filename
    os.system(cmd)
        
def determinePointOfSplit(options,chromosome_to_length):
    """
    Determines the location of split
    """
    chromosome_to_regions_with_N=findRegionsOfNs(options,chromosome_to_length)
    cmd="gt gff3 -retainids yes -addintrons yes "+options.gene_annotation+" > "+options.gene_annotation.split(".gff3")[0]+".introns_added.gff3"
    print(options.gene_annotation,options.gene_annotation.split(".gff3"))
    print(cmd)
    os.system(cmd)
    chromosome_to_split={}
    chromosome_to_regions_with_annotation={}
    fhr=open(options.gene_annotation.split(".gff3")[0]+".introns_added.gff3","r")
    for line in fhr:
        if "#" in line:continue
        if "gene" not in line:continue
        line=line.strip().split("\t")
        if "gene" in line[2]:
            if line[0] not in chromosome_to_regions_with_annotation:
                chromosome_to_regions_with_annotation[line[0]]=[]
            chromosome_to_regions_with_annotation[line[0]].append([int(line[3]),int(line[4])])
    fhr.close()
    print(chromosome_to_regions_with_annotation.keys())
    
    for chromosome in chromosome_to_regions_with_N:
        regions_with_N=chromosome_to_regions_with_N[chromosome]
        regions_with_annotation=chromosome_to_regions_with_annotation[chromosome]
        for eachregion_N in regions_with_N:
            start,end=eachregion_N
            if end-start+1<100:continue
            f=0
            for eachregion_annotation in regions_with_annotation:
                start_annotated,end_annotated=eachregion_annotation
                if start_annotated<math.pow(2,28):continue
                if not ((end_annotated<start and start_annotated<start) or (start_annotated>end and end_annotated>end)):
                    f=1
                    break
            if f==0:
                chromosome_to_split[chromosome]=int((start+end)/2)
                break
    #print(chromosome_to_split)
    return chromosome_to_split

def splitGenome(options,chromosome_to_split_location):
    """
    Splits the chromosomes 
    """
    fhw=open(options.out_prefix+"_split.fa","w")
    fhr=open(options.genome,"r")
    for line in fhr:
        if ">" in line:
            chromosome=line.strip().split()[0][1:]
            if chromosome not in chromosome_to_split_location:
                fhw.write(line+fhr.readline())
            else:
                split_location=chromosome_to_split_location[chromosome]
                sequence=fhr.readline().strip()
                sequence_A, sequence_B = sequence[:split_location-1], sequence[split_location-1:]
                print("Splitting up ",chromosome)
                fhw.write(">"+chromosome+"_A\n"+sequence_A+"\n")
                fhw.write(">"+chromosome+"_B\n"+sequence_B+"\n")
    fhr.close()
    fhw.close()

def splitGFF3(options,chromosome_to_split_location):
    """
    """
    print("Splitting up GFF3 file")
    fhw=open(options.out_prefix+"_split.gff3","w")
    fhr=open(options.gene_annotation,"r")
    for line in fhr:
        if "#" in line:
            if "sequence-region" not in line:
                fhw.write(line)
        else:
            chromosome=line.strip().split()[0]
            if chromosome not in chromosome_to_split_location:
                fhw.write(line)
            else:
                if int(line.strip().split("\t")[3]) < chromosome_to_split_location[line.strip().split("\t")[0]]:
                    line=line.strip().split("\t")
                    line[0]=line[0]+"_A"
                    fhw.write("\t".join(line)+"\n")
                else:
                    #print("Entered")
                    line=line.strip().split("\t")
                    line[3]=str(int(line[3])-(chromosome_to_split_location[line[0]]-1))
                    line[4]=str(int(line[4])-(chromosome_to_split_location[line[0]]-1))
                    line[0]=line[0]+"_B"
                    line="\t".join(line)+"\n"
                    fhw.write(line)
    fhr.close()
    fhw.close()
    cmd="gt gff3 -retainids yes -addids yes -tidy yes -addintrons yes -sort yes -fixregionboundaries yes "+options.out_prefix+"_split.gff3 > "+options.out_prefix+"_split.gff3.temp "
    os.system(cmd)
    cmd="mv "+options.out_prefix+"_split.gff3.temp "+options.out_prefix+"_split.gff3 "
    os.system(cmd)
    
def performSplit(options,chromosome_to_split_location):
    """
    Perform the actual split
    """
    splitGenome(options,chromosome_to_split_location)
    splitGFF3(options,chromosome_to_split_location)

def verifyForCorrectness(options,chromosome_to_split_location):
    """
    Looks into the gff3 file to ensure that the split region does not
    have any annotation
    """  
    #print(chromosome_to_split_location)
    for chromosome in chromosome_to_split_location:
        split_point=chromosome_to_split_location[chromosome]
        cmd="cat "+options.gene_annotation.split(".gff3")[0]+".introns_added.gff3 | "
        cmd+="awk '$1==\""+chromosome+"\"'|awk '$4<="+str(split_point)+" && $5>="+str(split_point)+"' > "
        cmd+=options.out_prefix+".temp"
        os.system(cmd)
        
        if len(open(options.out_prefix+".temp","r").read())!=0:
            return -1

def cleanUpFiles(options):
    """
    Remove all intermediate files
    """
    cmd="rm "
    cmd+=options.genome+".fai "
    #cmd+=options.out_prefix+"_regions_with_N.bed "
    cmd+=options.gene_annotation.split(".gff3")[0]+".introns_added.gff3 "
    #cmd+=options.gene_annotation.split(".gff3")[0]+".introns_added.bed "
    #cmd+=options.out_prefix+"_regions_with_N_covered_by_genes.bed "
    #cmd+=options.out_prefix+"_final_regions.bed "
    os.system(cmd)
    
def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    chromosome_to_length=findChromosomes(options)
    chromosome_to_split_location=determinePointOfSplit(options,chromosome_to_length)
    if verifyForCorrectness(options,chromosome_to_split_location)==-1:
        print("ERROR")
        sys.exit()
    performSplit(options,chromosome_to_split_location)
    #cleanUpFiles(options)

if __name__ == "__main__":
    main()