import sys
from pathlib import Path

# 获取当前文件的绝对路径
#current_file = Path(__file__).resolve()
# 获取 utils 目录的父目录（例如 pages 目录）
#parent_dir = current_file.parent.parent
# 将父目录添加到 Python 路径
#sys.path.append(str(parent_dir))
from .artfeat.ProteinDescriptors import toAAC,toDPC,toBLOSUM62,toASDC,toCTD,toAPAAC,toPAAC
from .artfeat.ToAAindex import toAAindex_emb
import pandas as pd
def get8artfeat(seq_pd):

    AAC_pd=toAAC(seq_pd)
    ASDC_pd=toASDC(seq_pd)
    DPC_pd=toDPC(seq_pd)
    PAAC_pd=toPAAC(seq_pd)
    AAindex_pd=toAAindex_emb(seq_pd)
    APAAC_pd=toAPAAC(seq_pd)
    CTD_pd=toCTD(seq_pd)
    BLOSUM62_pd=toBLOSUM62(seq_pd)
    features_pd=pd.concat([AAC_pd,ASDC_pd.iloc[:,2:],\
                           DPC_pd.iloc[:,2:],\
                           PAAC_pd.iloc[:,2:],\
                           AAindex_pd.iloc[:,2:],\
                           APAAC_pd.iloc[:,2:],\
                           CTD_pd.iloc[:,2:],\
                           BLOSUM62_pd.iloc[:,2:],\
                        
                           ],axis=1,ignore_index=True)
    #features_pd.to_csv("tempfeats.csv",index=False)
    #print(features_pd.shape)
    return features_pd
   
