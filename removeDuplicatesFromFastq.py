import sys

def readUniqueReadsFromFastqFile(filename):
    """
    Reads in the unique reads
    """
    fhr=open(filename,"r")
    data={}
    while True:
        read_id=fhr.readline().strip()
        if not read_id: break
        if read_id not in data:
            seq=fhr.readline().strip()
            useless=fhr.readline().strip()
            quality=fhr.readline().strip()
            data[read_id]={"seq":seq,"useless":useless,"quality":quality}
        else:
            seq=fhr.readline().strip()
            useless=fhr.readline().strip()
            quality=fhr.readline().strip()
    fhr.close()
    return data

def writeUniqueReadsToFastqFile(filename,unique_reads):
    """
    Write out the reads to file
    """
    fhw=open(filename,"w")
    for read_id in unique_reads:
        fhw.write(read_id+"\n"+unique_reads[read_id]["seq"]+"\n"+unique_reads[read_id]["useless"]+"\n"+unique_reads[read_id]["quality"]+"\n")
    fhw.close()

def main():
    if len(sys.argv)!=3:
        print("Usage removeDuplicates.py <name_of_fastq_file> <output>")
        sys.exit()
    unique_reads=readUniqueReadsFromFastqFile(sys.argv[1])
    writeUniqueReadsToFastqFile(sys.argv[2],unique_reads)
    
    
if __name__ == "__main__":
    main()