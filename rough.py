
"""for alpha in ["A","B","C","D","E","F","G","H"]:
    for num in range(1,13):
        if alpha=="H" and num>6:break
        cmd="nohup python /work/LAS/rpwise-lab/sagnik/scripts/retainBestReads.py "
        cmd+=" -i "+" /work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/trimmed/"+alpha+str(num)+".fastq "
        cmd+=" -p  /work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/trimmed/"+alpha+str(num)+" -u & "
        print(cmd)"""
        
filename="/work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/trimmed/all_combined_best.fastq"
fhr=open(filename,"r")
fhw=open(filename+".temp","w")
for line in fhr:
    fhw.write(line)

fhw.close()
fhr.close()