"""
Usage replaceGenomeSequence.py <name of the genome file> <chr2H:start-end> <sequence with which replacement needs to be made> <output_filename>
"""

import sys
genome_filename=sys.argv[1]
sequence_to_replaced=sys.argv[2]
replacement_seq=sys.argv[3]
output_filename=sys.argv[4]

chromosome=sequence_to_replaced.split(":")[0]
start=int(sequence_to_replaced.split(":")[-1].split("-")[0])
end=int(sequence_to_replaced.split(":")[-1].split("-")[1])
fhr=open(genome_filename,"r")
fhw=open(output_filename,"w")
while True:
    line=fhr.readline().strip()
    if not line : break
    if line[0]==">":
        #print(line.split()[0][1:],chromosome)
        if line.split()[0][1:]==chromosome:
            print("Making a change")
            fhw.write(line+"\n")
            seq=fhr.readline().strip()
            new_seq=seq[:start-1]+replacement_seq+seq[end-1:]
            fhw.write(new_seq+"\n")
        else:
            fhw.write(line+"\n")
            fhw.write(fhr.readline())
fhw.close()
fhr.close()