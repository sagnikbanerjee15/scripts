from collections import Counter

def readFasta(filename,type):
    data={}
    if type=="transcriptome":
        fhr=open(filename,"r")
        for line in fhr:
            if ">" in line:
                if "CDS=" in line:
                    seq=fhr.readline().strip()
                    data[line.strip()[1:].split()[0]]={"seq":seq,"cds":[int(line.strip().split()[-1].split("CDS=")[-1].split("-")[0]),
                                                                                           int(line.strip().split()[-1].split("CDS=")[-1].split("-")[-1])],
                                                                                    "5_prime_UTR_present":(0 if int(line.strip().split()[-1].split("CDS=")[-1].split("-")[0])==1 else 1),
                                                                                    "3_prime_UTR_present":(0 if int(line.strip().split()[-1].split("CDS=")[-1].split("-")[-1])==len(seq) else 1),
                                                                                    "5_prime_UTR_length":int(line.strip().split()[-1].split("CDS=")[-1].split("-")[0])-1,
                                                                                    "3_prime_UTR_length":len(seq)-int(line.strip().split()[-1].split("CDS=")[-1].split("-")[-1])
                    }
        fhr.close()
    else:
        fhr=open(filename,"r")
        for line in fhr:
            if ">" in line:
                seq=fhr.readline().strip()
                data[line.strip()[1:].split()[0]]={"seq":seq,"starts_with_M":(0 if seq[0]!="M" else 1)}
        fhr.close()
    return data

def readExpressionDataForGenes(filename):
    data={}
    fhr=open(filename,"r")
    for line in fhr:
        if "WT_0" in line:continue
        data[line.split(",")[0]]=[int(line.split(",")[15]),int(line.split(",")[16]),int(line.split(",")[17])]
    fhr.close()
    return data

def readTranscriptData(filename):
    transcriptome=readFasta(filename,"transcriptome")
    gene_expression_data=readExpressionDataForGenes("/work/LAS/rpwise-lab/sagnik/RNA-Seq_90_sample/data/salmon/90_samples_barley_bgh_genes.csv")
    proteome=readFasta("/work/LAS/rpwise-lab/sagnik/data/barley/proteome/Hordeum_vulgare.IBSC_v2.pep.all.fa", "protein")
    #print(len(transcriptome))
    genes_to_transcripts={}
    for transcript_id in transcriptome:
        if "HORVU" not in transcript_id: continue
        gene_id=transcript_id.split(".")[0]
        if gene_id not in genes_to_transcripts:
            genes_to_transcripts[gene_id]=[]
        genes_to_transcripts[gene_id].append(transcript_id)
    
    #Distribution of 5_prime_UTR Length
    """lengths=[]
    for gene_id in genes_to_transcripts:
        flag=0
        for transcript_id in genes_to_transcripts[gene_id]:
            lengths.append(transcriptome[transcript_id]["5_prime_UTR_length"]+transcriptome[transcript_id]["3_prime_UTR_length"])"""
    
    #print(dict(Counter(lengths)))
    
    """for key in sorted(dict(Counter(lengths)).keys()):
        print(key,dict(Counter(lengths))[key])"""
        
    starting_amino_acid=[]
    for gene_id in genes_to_transcripts:
        for transcript_id in genes_to_transcripts[gene_id]:
            starting_amino_acid.append(proteome[transcript_id]["seq"][0])
    
    for key in sorted(dict(Counter(starting_amino_acid)).keys()):
        print(key,dict(Counter(starting_amino_acid))[key])
    #Genes where all transcripts are missing 5 prime UTRs:
    """for gene_id in genes_to_transcripts:
        flag=0
        for transcript_id in genes_to_transcripts[gene_id]:
            if transcriptome[transcript_id]["5_prime_UTR_present"]==0:
                #print(transcript_id)
                print(transcript_id,gene_id,len(genes_to_transcripts[gene_id]),gene_expression_data[gene_id])
                flag=1"""
                #break
    """if flag==0:
            print("transcript:"+gene_id+".1",len(genes_to_transcripts[gene_id]),gene_expression_data[gene_id])"""
    
    #Genes where all transcripts are missing both UTRs:
    """for gene_id in genes_to_transcripts:
        flag=0
        for transcript_id in genes_to_transcripts[gene_id]:
            if transcriptome[transcript_id]["5_prime_UTR_present"]==1:
                flag=1
                print("transcript:"+transcript_id,gene_id,len(genes_to_transcripts[gene_id]),gene_expression_data[gene_id])
                #break"""
    """if flag==0:
            if sum(gene_expression_data[gene_id])>0:
                print("transcript:"+gene_id+".1",len(genes_to_transcripts[gene_id]),gene_expression_data[gene_id])"""
    
    
    """# Genes where no transcript begin with M
    for gene_id in genes_to_transcripts:
        flag=0
        for transcript_id in genes_to_transcripts[gene_id]:
            if proteome[transcript_id]["starts_with_M"]==1:
                flag=1
                break
        if flag==0:
            print("transcript:"+gene_id+".1",len(genes_to_transcripts[gene_id]),gene_expression_data[gene_id])"""
            
    """# Transcripts which do not start with M
    for gene_id in genes_to_transcripts:
        for transcript_id in genes_to_transcripts[gene_id]:
            if proteome[transcript_id]["starts_with_M"]==0:
                print(transcript_id,gene_id,len(genes_to_transcripts[gene_id]),gene_expression_data[gene_id])"""
    
    #num_of_transcripts=[len(genes_to_transcripts[gene_id]) for gene_id in genes_to_transcripts]
    """print(dict(Counter(num_of_transcripts)))
    for key in sorted(dict(Counter(num_of_transcripts)).keys()):
        print(key,dict(Counter(num_of_transcripts))[key])"""
    
    #print(",".join(list(map(str,num_of_transcripts))))
    #print(num_of_transcripts)
    
            

def main():
    readTranscriptData("/work/LAS/rpwise-lab/sagnik/data/barley/transcriptome/transcripts_gffread.fasta")

if __name__ == "__main__":
    main()