import pandas as pd
import numpy as np 
from collections import Counter
import os, sys, re
import itertools
import math
import random
import pickle
from rich.progress import track
import time 
        

def Count( seq1, seq2):
    sum = 0
    for aa in seq1:
        sum = sum + seq2.count(aa)
    return sum

    
def Count1(aaSet, sequence):
    number = 0
    for aa in sequence:
        if aa in aaSet:
            number = number + 1
    cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
    cutoffNums = [i if i >= 1 else 1 for i in cutoffNums]

    code = []
    for cutoff in cutoffNums:
        myCount = 0
        for i in range(len(sequence)):
            if sequence[i] in aaSet:
                myCount += 1
                if myCount == cutoff:
                    code.append((i + 1) / len(sequence) * 100)
                    break
        if myCount == 0:
            code.append(0)
    return code
        

def Rvalue(aa1, aa2, AADict, Matrix):
        return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)

def generatePropertyPairs(myPropertyName):
    pairs = []
    for i in range(len(myPropertyName)):
        for j in range(i + 1, len(myPropertyName)):
            pairs.append([myPropertyName[i], myPropertyName[j]])
            pairs.append([myPropertyName[j], myPropertyName[i]])
    return pairs
        
def toAAC(seq_pd):
     
    encoding_array = np.array([])
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    header=[]
    for aa in AA:
        header.append("Freq"+aa)
    encodings = []
    
    
    SEQ_=seq_pd["Sequence"]
 
    PID_=seq_pd["ID"]


    
    for sequence in track(SEQ_,"Computing: "):
        count = Counter(sequence)
        for key in count:
            count[key] = count[key] / len(sequence)
        code = []
        for aa in AA:
            code.append(count[aa])
        
        #print(code)
        encodings.append(code)
        
    encoding_array = np.array(encodings, dtype=str)
    AAC_encoding=pd.DataFrame(encoding_array)
    AAC_encoding.columns=header
    
    AAC_encoding=pd.concat([seq_pd,AAC_encoding],axis=1,ignore_index=True)
    
    return AAC_encoding
    
def toDPC(seq_pd):
    encoding_array = np.array([])

    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    diPeptides = ["Freq_"+aa1 + aa2 for aa1 in AA for aa2 in AA]
   
    SEQ_=seq_pd["Sequence"]
     
   
    PID_=seq_pd["ID"]
    #print(header)

    DPC_list=[aa1 + aa2 for aa1 in AA for aa2 in AA]
    
    #print(DPC_list)
    
    

    
    for sequence in track(SEQ_,"computing:"):
        #name, sequence, Class = i[0], re.sub("-","",i[1]), str(i[2])
        #code = [name, Class]
        tmpCode = [0] * 400
        
        dpc_seq=[sequence[k]+sequence[k+1] for k in range(len(sequence)-1)]
        
        #print(dpc_seq)
        
        count = Counter(dpc_seq)
        
        for key in count:
            count[key] = count[key] / (len(sequence)-1)
            
        code =[]
        for idpc in DPC_list:
            code.append(count[idpc])
        
        #print(code)
        encodings.append(code)
        
        
        
         
     

    encoding_array = np.array(encodings, dtype=str)
    
    DPC_encoding=pd.DataFrame(encodings)
    
    DPC_encoding.columns=diPeptides

    DPC_encoding=pd.concat([seq_pd,DPC_encoding],axis=1,ignore_index=True)
    
  
    return DPC_encoding
    
def toCTDC(seq_pd):
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity': 'LIFWCMVY',
        'polarizability': 'GASDT',
        'charge': 'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess': 'ALFCGIVW'
        }
            
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity': 'PATGS',
        'polarizability': 'CPNVEQIL',
        'charge': 'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess': 'RKQEND'
        }
    
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity': 'HQRKNED',
        'polarizability': 'KMHFRYW',
        'charge': 'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess': 'MSPTHY'
    }

    groups = [group1, group2, group3]
    property = (
                'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
                'hydrophobicity_PONP930101',
                'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
                'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []
    header = []
    
    for p in property:
        for g in range(1, len(groups) + 1):
            header.append(p + '.G' + str(g)+"_CTD.C")
            
            
    
     
    SEQ_=seq_pd["Sequence"]
     
    
    PID_=seq_pd["ID"]
    
    for sequence in track(SEQ_,"computing...."):
        code = []
        for p in property:
        
            c1 = Count(group1[p], sequence) / len(sequence)
            c2 = Count(group2[p], sequence) / len(sequence)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]
        encodings.append(code)

    encoding_array = np.array(encodings, dtype=str)
    
    CTDC_encoding=pd.DataFrame(encodings)
    
    CTDC_encoding.columns=header
    
    CTDC_encoding=pd.concat([seq_pd,CTDC_encoding],axis=1,ignore_index=True)
    CTDC_encoding.index=PID_
     
    
    
    return CTDC_encoding
    
def toCTDT(seq_pd):
    group1 = {
                'hydrophobicity_PRAM900101': 'RKEDQN',
                'hydrophobicity_ARGP820101': 'QSTNGDE',
                'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
                'hydrophobicity_PONP930101': 'KPDESNQT',
                'hydrophobicity_CASG920101': 'KDEQPSRNTG',
                'hydrophobicity_ENGD860101': 'RDKENQHYP',
                'hydrophobicity_FASG890101': 'KERSQD',
                'normwaalsvolume': 'GASTPDC',
                'polarity': 'LIFWCMVY',
                'polarizability': 'GASDT',
                'charge': 'KR',
                'secondarystruct': 'EALMQKRH',
                'solventaccess': 'ALFCGIVW'
            }
    group2 = {
                'hydrophobicity_PRAM900101': 'GASTPHY',
                'hydrophobicity_ARGP820101': 'RAHCKMV',
                'hydrophobicity_ZIMJ680101': 'HMCKV',
                'hydrophobicity_PONP930101': 'GRHA',
                'hydrophobicity_CASG920101': 'AHYMLV',
                'hydrophobicity_ENGD860101': 'SGTAW',
                'hydrophobicity_FASG890101': 'NTPG',
                'normwaalsvolume': 'NVEQIL',
                'polarity': 'PATGS',
                'polarizability': 'CPNVEQIL',
                'charge': 'ANCQGHILMFPSTWYV',
                'secondarystruct': 'VIYCWFT',
                'solventaccess': 'RKQEND'
            }
    group3 = {
                'hydrophobicity_PRAM900101': 'CLVIMFW',
                'hydrophobicity_ARGP820101': 'LYPFIW',
                'hydrophobicity_ZIMJ680101': 'LPFYI',
                'hydrophobicity_PONP930101': 'YMFWLCVI',
                'hydrophobicity_CASG920101': 'FIWC',
                'hydrophobicity_ENGD860101': 'CVLIMF',
                'hydrophobicity_FASG890101': 'AYHWVMFLIC',
                'normwaalsvolume': 'MHKFRYW',
                'polarity': 'HQRKNED',
                'polarizability': 'KMHFRYW',
                'charge': 'DE',
                'secondarystruct': 'GNPSD',
                'solventaccess': 'MSPTHY'
            }

    groups = [group1, group2, group3]
    property = (
                'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
                'hydrophobicity_PONP930101',
                'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
                'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []
    header = []
    
    for p in property:
        for tr in ('Tr1221', 'Tr1331', 'Tr2332'):
            header.append(p + '.' + tr)
    
     
    SEQ_=seq_pd["Sequence"]
     
     
    PID_=seq_pd["ID"]

    for sequence in track(SEQ_,"computing..."):
       
        code = []
        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
        
        for p in property:
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 = c1221 + 1
                    continue
                    
                if (pair[0] in group1[p] and pair[1] in group3[p]) or ( pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 = c1331 + 1
                    continue
                    
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 = c2332 + 1
                
            code = code + [c1221 /(1e-6+ len(aaPair)), c1331 / (1e-6+len(aaPair)), c2332 /(1e-6+ len(aaPair))]
        encodings.append(code)

    #encoding_array = np.array(encodings, dtype=str)
    
    CTDT_encoding=pd.DataFrame(encodings)
    
    CTDT_encoding.columns=header
    
    CTDT_encoding=pd.concat([seq_pd,CTDT_encoding],axis=1,ignore_index=True)
        
     
    return CTDT_encoding

def toCTDD(seq_pd):

    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity': 'LIFWCMVY',
        'polarizability': 'GASDT',
        'charge': 'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess': 'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity': 'PATGS',
        'polarizability': 'CPNVEQIL',
        'charge': 'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess': 'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity': 'HQRKNED',
        'polarizability': 'KMHFRYW',
        'charge': 'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess': 'MSPTHY'
    }

    groups = [group1, group2, group3]
    property = (
        'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101',
        'hydrophobicity_PONP930101',
        'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
        'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    encodings = []
    header = []
    for p in property:
        for g in ('1', '2', '3'):
            for d in ['0', '25', '50', '75', '100']:
                header.append(p + '.' + g + '.residue' + d)
    
    
     
    SEQ_=seq_pd["Sequence"]
     
     
    PID_=seq_pd["ID"]

    for sequence in track(SEQ_,"computing..."):

        #name, sequence, Class = i[0], re.sub('-', '', i[1]), str(i[2])
        code = []
        for p in property:
            code = code +  Count1(group1[p], sequence) +  Count1(group2[p], sequence) +  Count1(
                group3[p], sequence)
                
        encodings.append(code)

    CTDD_encoding=pd.DataFrame(encodings)
    
    CTDD_encoding.columns=header
    
    CTDD_encoding=pd.concat([seq_pd,CTDD_encoding],axis=1,ignore_index=True)
    
    return CTDD_encoding
    
def toCTD(seq_pd):
    
    ctdc=toCTDC(seq_pd)
    print(ctdc.shape)
    ctdt=toCTDT(seq_pd)
    print(ctdt.shape)
    ctdd=toCTDD(seq_pd)
    print(ctdd.shape)

    # 先重置索引（关键！）
    ctdc = ctdc.reset_index(drop=True)
    ctdt = ctdt.reset_index(drop=True)
    ctdd = ctdd.reset_index(drop=True)

# 再横向合并（列合并）
    ctd = pd.concat([ctdc, ctdt.iloc[:,2:], ctdd.iloc[:,2:]], axis=1)
    
    #ctd=pd.concat([ctdc,ctdt.iloc[:,2:],ctdd.iloc[:,2:]],axis=1,ignore_index=True)
     
   
    return ctd
    

def toPAAC(seq_pd,lambdaValue = 2,w = 0.05):
 
    lambdaValue = lambdaValue
    w = w
    dataFile = os.path.join(os.path.dirname(__file__), 'util_data', 'PAAC.txt')
    
    with open(dataFile) as f:
        records = f.readlines()
        
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    
    for i in range(len(AA)):
        AADict[AA[i]] = i
        
    AAProperty = []
    AAPropertyNames = []
    
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])

    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])
    
    encodings = []
    header = []
    
    for aa in AA:
        header.append('Xc1.' + aa)
    for n in range(1, lambdaValue + 1):
        header.append('Xc2.lambda' + str(n))
     
     
    SEQ_=seq_pd["Sequence"]
     
    
    PID_=seq_pd["ID"]
    
    for sequence in track(SEQ_,"PAAC Computing..."):
        #name, sequence, Class = i[0], re.sub('-', '', i[1]), str(i[2])
        code = []
        
        theta = []
        for n in range(1, lambdaValue + 1):
            if (len(sequence) - n)==0:

                theta.append(
                    sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in
                        range(len(sequence) - n)]) / (1+
                            len(sequence) - n))
            else:
                theta.append(
                    sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in
                        range(len(sequence) - n)]) / (
                            len(sequence) - n))

        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)
        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
        encodings.append(code)
        
    
    PAAC_encoding=pd.DataFrame(encodings)
    
    PAAC_encoding.columns=header
    PAAC_encoding=pd.concat([seq_pd,PAAC_encoding],axis=1,ignore_index=True)
     
    
    return PAAC_encoding  
    

def toAPAAC(seq_pd,lambdaValue = 2,weight = 0.05):
 
    lambdaValue =  lambdaValue 
    
    w = weight
    dataFile= os.path.join(os.path.dirname(__file__), 'util_data', 'PAAC.txt')
    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records) - 1):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])

    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])

    encodings = []
    header = []
    for i in AA:
        header.append('Pc1.' + i)
    for j in range(1, lambdaValue + 1):
        for i in AAPropertyNames:
            header.append('Pc2.' + i + '.' + str(j))
    
    
     
    SEQ_=seq_pd["Sequence"]
     
     
    PID_=seq_pd["ID"]
    
    for sequence in  track(SEQ_,"APAAC computing..."):
        #name, sequence, Class = i[0], re.sub('-', '', i[1]), str(i[2])
        code = []
        theta = []
        for n in range(1, lambdaValue + 1):
            for j in range(len(AAProperty1)):
                if (len(sequence) - n)==0:
                    theta.append(
                    sum([AAProperty1[j][AADict[sequence[k]]] * AAProperty1[j][AADict[sequence[k + n]]] for k in
                        range(len(sequence) - n)]) / (1+len(sequence) - n))

                else:            

                    theta.append( sum([AAProperty1[j][AADict[sequence[k]]] * AAProperty1[j][AADict[sequence[k + n]]] for k in
                        range(len(sequence) - n)]) / (len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)

        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [w * value / (1 + w * sum(theta)) for value in theta]
        encodings.append(code)
    
    APAAC_encoding=pd.DataFrame(encodings)
    
    APAAC_encoding.columns=header
    
    APAAC_encoding=pd.concat([seq_pd,APAAC_encoding],axis=1,ignore_index=True)
    
    return APAAC_encoding


     
def toASDC(seq_pd):
     
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    aaPairs = []
    for aa1 in AA:
        for aa2 in AA:
            aaPairs.append(aa1 + aa2)

    header = []
    header += ["ASDC_"+aa1 + aa2 for aa1 in AA for aa2 in AA]
    
     
    SEQ_=seq_pd["Sequence"]
     
     
    PID_=seq_pd["ID"]

    for sequence in track(SEQ_,"ASDC computing..."):
        #name, sequence, Class = i[0], re.sub('-', '', i[1]), str(i[2])
        code = []
        sum = 0
        pair_dict = {}
        for pair in aaPairs:
            pair_dict[pair] = 0
        for j in range(len(sequence)):
            for k in range(j + 1, len(sequence)):
                if sequence[j] in AA and sequence[k] in AA:
                    pair_dict[sequence[j] + sequence[k]] += 1
                    sum += 1
        for pair in aaPairs:
            code.append(pair_dict[pair] /(1e-6+ sum))
        encodings.append(code)
    
    ASDC_encoding=pd.DataFrame(encodings)
    
    ASDC_encoding.columns=header
    ASDC_encoding=pd.concat([seq_pd,ASDC_encoding],axis=1,ignore_index=True)
     
    
    return ASDC_encoding
     
     
def toBLOSUM62(seq_pd):
 
    
    blosum62 = {
        'A': [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],  # A
        'R': [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
        'N': [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],  # N
        'D': [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],  # D
        'C': [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
        'Q': [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],  # Q
        'E': [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],  # E
        'G': [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],  # G
        'H': [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],  # H
        'I': [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],  # I
        'L': [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],  # L
        'K': [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],  # K
        'M': [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],  # M
        'F': [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],  # F
        'P': [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],  # P
        'S': [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],  # S
        'T': [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],  # T
        'W': [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],  # W
        'Y': [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],  # Y
        'V': [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],  # V
        '-': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # -
    }
    encodings = []
    header = []
    
    AA='ACDEFGHIKLMNPQRSTVWY-'
    
    
   
    
    for i in range(21*20):
        header.append('blosum62.F' + str(i+1))
    
    
    SEQ_=seq_pd["Sequence"]
     
    
    PID_=seq_pd["ID"]
   

    for sequence in track(SEQ_," blosum62 computing..."):
        #name, sequence, Class = i[0], i[1], i[2]
        code = []
        count=Counter(sequence)
        #print(count)
       # print(count.keys())
        for aa in AA:
            if aa in count.keys():
                code = code + [i*count[aa] for i in blosum62[aa]]
            else:
                code = code +  blosum62['-']
            
        encodings.append(code)
    
    
    blosum62_encoding=pd.DataFrame(encodings)
    #print(blosum62_encoding)
    

    blosum62_encoding.columns=header
    blosum62_encoding=pd.concat([seq_pd,blosum62_encoding],axis=1,ignore_index=True)
     

    return blosum62_encoding

     
    
