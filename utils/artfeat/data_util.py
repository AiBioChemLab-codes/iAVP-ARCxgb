import pandas as pd
from Bio import SeqIO
def read_fasta(inFasta):
    record_ids=[]
    record_seqs=[]
    record_targets=[]
    for seq_record in SeqIO.parse(inFasta, "fasta"):
        #print(seq_record.id)
        tem=seq_record.id.split("|")
        #print()
        idss=""
        for ss in range(len(tem)-1):
          idss=idss+str(tem[ss])
        #print(len(tem),idss)
          
        record_ids.append(idss)
        record_targets.append(str(tem[-1]))
        #print(seq_record.seq)
        record_seqs.append(str(seq_record.seq).upper())
    outpd=pd.DataFrame()
    outpd["PID"]=record_ids
    outpd["Target"]=record_targets
    outpd["Sequence"]=record_seqs
    return outpd

def assert_filetype(filename):
    fn=filename.split(".")[-1]
    if str(fn).lower()=="csv":
      return "csv"
    elif str(fn).lower() in ["fa","fasta",'txt']:
      return "fasta"
    else:
      print("不是csv或fasta格式!")
      exit()

def read_seq_file(infile):
    file_type= assert_filetype(infile)
    if file_type=="csv":
        data=pd.read_csv(infile,header=0)
        
    elif file_type=="fasta":
        data=read_fasta(infile)
    else:
        print("请按照要求输入指定格式的文件csv格式或fasta格式!")
        exit()
    return data

