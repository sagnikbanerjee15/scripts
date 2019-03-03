from argparse import ArgumentParser
import sys

def extractReads(sam_filename,read_ids_filename,output_filename):
    """
    """
    fhr=open(read_ids_filename,"r")
    read_ids=set([line.strip().split()[0] for line in fhr])
    fhr.close()
    fhr=open(sam_filename)
    fhw=open(output_filename,"w")
    prev_read_id=""
    found=0
    for line in fhr:
        if line[0]=="@":
            fhw.write(line)
            continue
        read_id=line.split()[0]
        if found==1 and prev_read_id==read_id:
            fhw.write(line)
        elif found==0 and prev_read_id==read_id:
            continue
        elif prev_read_id!=read_id:
            if read_id in read_ids:
                fhw.write(line)
                prev_read_id=read_id
                found=1
    fhw.close()
            

def main():
    if len(sys.argv)!=3:
        print("Please enter the sam file and a newline separated read ids file.\n For example python extractReadsFromNameSortedSamFile.py xyz.sam abc.txt")
    sam_filename=sys.argv[1]
    read_ids_filename=sys.argv[2]
    output_filename=sam_filename.split(".sam")[0]+"_subset.sam"
    extractReads(sam_filename,read_ids_filename,output_filename)

if __name__ == "__main__":
    main()