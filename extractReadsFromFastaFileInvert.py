
import sys
import os


def extractReadsInvert(fastafilename,readfilenames):
    """
    """
    all_reads=[]
    for filename in readfilenames:
        all_reads.extend([line.strip() for line in open(filename,"r")])
    all_reads=set(all_reads)
    #print(list(all_reads)[:10])
    fhr=open(fastafilename,"r")
    while True:
        id=fhr.readline()
        #print(id[1:].strip().split()[0])
        if not line:break
        if id[1:].strip().split()[0] in all_reads:
            fhr.readline()
        else:
            print(id.strip())
            print(fhr.readline().strip())
    fhr.close()
    
def main():
    if len(sys.argv)<1:
        print("Usage: python extractReadsFromFastqFileInvert.py <fastafilename> Output will be on stdout")
    fastafilename=sys.argv[1]
    readfilenames=[sys.argv[2]]
    extractReadsInvert(fastafilename,readfilenames)

if __name__ == "__main__":
    main()