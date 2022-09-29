
import sys
import os


def extractReads(fastqfilename,readfilenames):
    """
    """
    all_reads=[]
    for filename in readfilenames:
        all_reads.extend([line.strip() for line in open(filename,"r")])
    all_reads=set(all_reads)
    #print(list(all_reads)[:10])
    fhr=open(fastqfilename,"r")
    while True:
        id=fhr.readline()
        #print(id[1:].strip().split()[0])
        if not line:break
        if id[1:].strip().split()[0] not in all_reads:
            fhr.readline()
            fhr.readline()
            fhr.readline()
        else:
            print(id.strip())
            print(fhr.readline().strip())
            print(fhr.readline().strip())
            print(fhr.readline().strip())
    fhr.close()
    
def main():
    if len(sys.argv)<1:
        print("Usage: python extractReadsFromFastqFile.py <fastqfilename> <readfilename1> <readfilename2> ... Output will be on stdout")
    fastqfilename=sys.argv[1]
    readfilenames=sys.argv[2:]
    extractReads(fastqfilename,readfilenames)

if __name__ == "__main__":
    main()