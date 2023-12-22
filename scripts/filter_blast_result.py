import sys
import pandas as pd
from Bio.Seq import Seq

def filter_by_identity_coverage(blast_file, thre_identity, thre_coverage,output_file):    
    #read file
    blast_result = pd.read_table(blast_file,sep="\t",header=None, 
                                 names= ["query", "subject", "identity","length","mismatch","gapopen", "qstart","qend","sstart","send","evalue","bitscore","gaps","qcovs","qcovhsp","sstrand","sseq"])
    if blast_result.empty:
        print("Empty blast result!")
        blast_result.to_csv(output_file,sep="\t",header=True, index=False)
    else:
        #create VFid
        blast_result["query_short"] = blast_result["query"].str.replace('\|.*\|[0-9]*-[0-9]*','',regex=True)
        #filter by identity and coverge
        blast_result_1 = blast_result[(blast_result['identity'] > float(thre_identity)) & (blast_result['qcovs'] > float(thre_coverage))]
        if blast_result_1.empty:
            print("No hits fullfills!")
            blast_result_1.to_csv(output_file,sep="\t",header=True, index=False)
        else:
            #select the best hit for subject region, by coverage, identity and bitscore
            blast_result_2 = blast_result_1.sort_values(['qcovs','identity','bitscore'],ascending=False).groupby(['query_short','subject','sstart','send']).first().reset_index()
            #evaluate pathogen alleles
            blast_result_2["type"] = blast_result_2.apply(lambda row : "Perfect" if (row["identity"]==100 & row["qcovs"]==100) else "Strict", axis=1)
            #transform seq in the minus strand
            blast_result_2["sseq_new"] = blast_result_2.apply(lambda row : row['sseq'].replace("-",'')  if (row["sstart"] < row["send"]) else ''.join(Seq(row['sseq'].replace("-",'')).reverse_complement()), axis=1)
            #transform pos in the minus strand
            blast_result_2["sstart_new"] = blast_result_2.apply(lambda row: row["sstart"] if (row["sstart"] < row["send"]) else row["send"], axis=1)
            blast_result_2["ssend_new"] = blast_result_2.apply(lambda row: row["send"] if (row["sstart"] < row["send"]) else row["sstart"], axis=1)
            #blast_result_2["sstrand"] = blast_result_2.apply(lambda row: "plus" if (row["sstart"] < row["send"]) else "minus", axis=1)
            blast_result_2.to_csv(output_file,sep="\t",header=True, index=False)


filter_by_identity_coverage(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
