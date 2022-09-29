import os
import sys

def readGTFFile(filename):
    """
    """
    gtf_info={"+":{},"-":{}}
    fhr=open(filename,"r")
    for line in fhr:
        line=line.strip().split("\t")
        if line[2]!="exon":continue
        gene_id=line[-1].split(";")[1].split("gene_id")[-1].strip().strip("\"")
        strand=line[6]
        start=int(line[3])
        end=int(line[4])
        chromosome=line[0]
        if chromosome not in gtf_info[strand]:
            gtf_info[strand][chromosome]={}
        if gene_id not in gtf_info[strand][chromosome]:
            gtf_info[strand][chromosome][gene_id]={"start":start,"end":end}
        else:
            if start < gtf_info[strand][chromosome][gene_id]["start"]:
                gtf_info[strand][chromosome][gene_id]["start"]=start
            if end > gtf_info[strand][chromosome][gene_id]["end"]:
                gtf_info[strand][chromosome][gene_id]["end"]=end
    fhr.close()
    return gtf_info

def findCloselyPlacedGenes(gene_info):
    """
    """
    

def main():
    if len(sys.argv)<2:
        print("python find_closely_spaced_genes.py <GTF file>. Closely spaced genes will be written to stdout. \n Run python /work/LAS/rpwise-lab/sagnik/scripts/find_closely_placed_genes.py /work/LAS/rpwise-lab/sagnik/data/arath/transcriptome/Arabidopsis_thaliana.TAIR10.42.gtf")
        sys.exit()
    gene_info=readGTFFile(sys.argv[1])
    """for strand in gene_info:
        for chromosome in gene_info[strand]:
            for gene in gene_info[strand][chromosome]:
                print(gene,chromosome,strand,gene_info[strand][chromosome][gene]["start"],gene_info[strand][chromosome][gene]["end"])"""
    findCloselyPlacedGenes(gene_info)

if __name__ == "__main__":
    main()