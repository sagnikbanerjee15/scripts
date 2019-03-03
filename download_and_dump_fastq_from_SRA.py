"""
Parallely downloads and converts NCBI .sra files to FASTQ
"""

import argparse
import os
import sys
from pathlib import Path
import multiprocessing

def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python download_and_dump_fastq_from_SRA.py --help for more details
    """
    parser = argparse.ArgumentParser(prog="download_and_dump_fastq_from_SRA.py",description="Parallel download of fastq data from NCBI.")
    parser.add_argument("--sra","-s",help="Please enter the name of the file which has all the SRAs listed one per line",required=True)
    parser.add_argument("--output","-o",help="Please enter the name of the output directory.",required=True)
    parser.add_argument("--cpu","-n",help="Enter the number of CPUs to be used.",default=1)
    return parser.parse_args()

def readSRAfilesToBeDownloaded(filename):
    """
    Reads and returns a list of the SRA ids to be downloaded
    """
    return list(set([name.strip() for name in open(filename,"r").read().split("\n")]))

def downloadSRAFile(sra,default_path_to_download,output_directory):
    os.system("prefetch "+sra)
    os.system("mv "+default_path_to_download+"/"+sra+".sra "+output_directory+"/")
    cmd="fastq-dump -X 1 -Z  --split-spot "+output_directory+"/"+sra+".sra|wc -l > "+output_directory+"/"+sra+".temp"
    os.system(cmd)
    if int(open(output_directory+"/"+sra+".temp").read())==4:
        pair="single"
    else:
        pair="paired"
    cmd="fastq-dump --defline-seq '@$sn[_$rn]/$ri' --outdir "+output_directory+" --split-files "+output_directory+"/"+sra+".sra"
    os.system(cmd)
    if pair=="single":
        os.system("mv "+output_directory+"/"+sra+"_1.fastq "+output_directory+"/"+sra+".fastq ")
    os.system("rm "+output_directory+"/"+sra+".sra "+output_directory+"/"+sra+".temp")
    
def downloadSRAFilesAndConvertToFastq(SRAs,default_path_to_download,n,output_directory):
    """
    Downloads the sra files and converts to fastq
    """
    cmd="mkdir "+output_directory
    os.system(cmd)
    jobs=[]
    for sra in SRAs:
        if os.path.exists(output_directory+"/"+sra+".fastq")==True or (os.path.exists(output_directory+"/"+sra+"_1.fastq")==True and os.path.exists(output_directory+"/"+sra+"_2.fastq")==True):
            if os.path.exists(output_directory+"/"+sra+"_1.fastq")==True and os.path.exists(output_directory+"/"+sra+"_2.fastq")==False:
                os.system("mv "+output_directory+"/"+sra+"_1.fastq "+output_directory+"/"+sra+".fastq")
            continue
        p = multiprocessing.Process(target=downloadSRAFile, args=(sra,default_path_to_download,output_directory))
        jobs.append(p)
    for p in jobs:
        while [p.is_alive() for p in jobs].count(True) > n: continue
        p.start()
    for p in jobs:p.join()
    
def verifyOutput(output_directory,SRAs):
    """
    Verify the downloads
    """    
    for sra in SRAs:
        if os.path.exists(output_directory+"/"+sra+".fastq")==True:
            continue
        elif os.path.exists(output_directory+"/"+sra+"_1.fastq")==True and os.path.exists(output_directory+"/"+sra+"_2.fastq")==True:
            continue
        print(sra,"was not downloaded. Please try manually.")
        
    
def main():
    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")
    options=parseCommandLineArguments()
    SRAs=readSRAfilesToBeDownloaded(options.sra)
    #home = str(Path.home())
    default_path_to_download="/home/sagnik/ncbi/public/sra/"
    downloadSRAFilesAndConvertToFastq(SRAs,default_path_to_download,int(options.cpu),options.output)
    verifyOutput(options.output,SRAs)

if __name__ == "__main__":
    main()
