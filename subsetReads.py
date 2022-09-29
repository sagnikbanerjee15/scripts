
import sys

def subsetReads(inputfilename,outputfilename,readfilename):
    all_reads=[line.strip() for line in open(readfilename)]
    all_reads=set(all_reads)
    print(all_reads[:5])
    sys.stdout.flush()
    fhr=open(inputfilename,"r")
    fhw=open(outputfilename,"w")
    sequences_read=0
    sequences_written=0
    for line_num,line in enumerate(fhr):
        if line_num%4==0:
            sequences_read+=1
            id=line.strip()
            if id.split()[0] in all_reads:
                print(id.split()[0])
                sys.stdout.flush()
                fhw.write(id+"\n"+fhr.readline()+fhr.readline()+fhr.readline())
                sequences_written+=1
                if sequences_written%1000==0:
                    print(sequences_read,sequences_written)
                    sys.stdout.flush()
    fhr.close()
    fhw.close()
        
subsetReads("/work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/trimmed/combined_best.fastq", 
            "/work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/trimmed/combined_best_blumeria.fastq", 
            "/work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/trimmed/STAR_alignments/blumeria_mapped_reads_and_unmapped_to_both")