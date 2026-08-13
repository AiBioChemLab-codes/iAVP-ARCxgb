import pandas as pd
import numpy as np 
from collections import Counter
import os, sys, re
import itertools
import math
import random
import pickle
import re
import platform
 


def read_fasta(iFasta):
    if not os.path.exists(iFasta):
        msg = 'Error: file %s does not exist.' % iFasta
        print(msg)
        return []
    with open(iFasta) as f:
        records = f.read()
        records = records.split('>')[1:]
        fasta_sequences = []
        
        for fasta in records:
            array = fasta.split('\n')
            header, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTUVWY-]', '-', ''.join(array[1:]).upper())
            header_array = header.split('|')
            name = header_array[0]
            label = header_array[1] if len(header_array) >= 2 else '0'
            #label_train = header_array[2] if len(header_array) >= 3 else 'training'
            fasta_sequences.append([name, sequence, label])
        msg="read_fasta DONE!"
        print(msg)
        return fasta_sequences 
        


    

def toNCP(inFasta,outCSV): #RNA和DNA都适用
    
    chemical_property = {
        'A': [1, 1, 1],
        'C': [0, 1, 0],
        'G': [1, 0, 0],
        'T': [0, 0, 1],
        'U': [0, 0, 1],
        '-': [0, 0, 0],
    }
     
    encodings = []
    header = ['PID', 'Label']
    
    fasta_list=read_fasta(inFasta)
    for i in range(1, len(fasta_list[0][1]) * 3 + 1):
        header.append('NCP.F' + str(i))
    

    for i in fasta_list:
        name, sequence, label = i[0], i[1], str(i[2])
        code = [name, label]
        for aa in sequence:
            code = code + chemical_property.get(aa, [0, 0, 0])
        encodings.append(code)
    
    encoding_array = np.array(encodings, dtype=str)
    NCP_encoding=pd.DataFrame(encoding_array)
    NCP_encoding.columns=header
    NCP_encoding.to_csv(outCSV,index=False)
    print("NCP is DONE. NCP_encoding=",NCP_encoding.shape)
    return NCP_encoding
     
def toBINARY(inFasta,outCSV,A_type="DNA"):#RNA和DNA都适用
    
    if A_type=="RNA":
        AA="ACGU"
    elif A_type=="DNA":
        AA = 'ACGT'
    encodings = []
    header = ['PID', 'Label']
    fasta_list=read_fasta(inFasta)
    
    for i in range(1, len(fasta_list[0][1]) * 4 + 1):
        header.append('BINARY.F' + str(i))
    
    
    
    for i in fasta_list:
        name, sequence, label = i[0], i[1], str(i[2])
        code = [name, label]
        for aa in sequence:
            if aa == '-':
                code = code + [0, 0, 0, 0]
                continue
            for aa1 in AA:
                tag = 1 if aa == aa1 else 0
                code.append(tag)
        encodings.append(code)
    
    encoding_array = np.array(encodings, dtype=str)
    BINARY_encoding=pd.DataFrame(encoding_array)
    BINARY_encoding.columns=header
    BINARY_encoding.to_csv(outCSV,index=False)
    print("BINARY is DONE. BINARY_encoding=",BINARY_encoding.shape)
    return BINARY_encoding
    

def toCKSNAP(inFasta,outCSV,gap=3,A_type="DNA"):
 
    if A_type=="RNA":
        AA="ACGU"
    elif A_type=="DNA":
        AA = 'ACGT'
        
    
    encodings = []
    aaPairs = []
    for aa1 in AA:
        for aa2 in AA:
            aaPairs.append(aa1 + aa2)

    header = ['PID', 'label']
    
    for g in range(gap + 1):
        for aa in aaPairs:
            header.append(aa + '.gap' + str(g))
    
    fasta_list=read_fasta(inFasta)

    for i in  fasta_list:
        name, sequence, label = i[0], i[1], str(i[2])
        code = [name, label]
        for g in range(gap + 1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(sequence)):
                index2 = index1 + g + 1
                if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[
                    index2] in AA:
                    myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                code.append(myDict[pair] / sum)
        encodings.append(code)
    
    encoding_array = np.array(encodings, dtype=str)
    CKSNAP_encoding=pd.DataFrame(encoding_array)
    CKSNAP_encoding.columns=header
    CKSNAP_encoding.to_csv(outCSV,index=False)
    print("CKSNAP is DONE. CKSNAP_encoding=",CKSNAP_encoding.shape)
    return CKSNAP_encoding

    
def toEIIP(inFasta,outCSV):
 
    fasta_list=read_fasta(inFasta)
    
    
    EIIP_dict = {
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'T': 0.1335,
        'U': 0.1335,
        '-': 0,
    }
    encodings = []
    header = ['PID', 'Label']
    for i in range(1, len(fasta_list[0][1]) + 1):
        header.append('F' + str(i))
     

    for i in fasta_list:
        name, sequence, label = i[0], i[1], str(i[2])
        code = [name, label]
        for aa in sequence:
            code.append(EIIP_dict.get(aa, 0))
        encodings.append(code)
        
    encoding_array = np.array(encodings, dtype=str)
    EIIP_encoding=pd.DataFrame(encoding_array)
    EIIP_encoding.columns=header
    EIIP_encoding.to_csv(outCSV,index=False)
    print("EIIP is DONE. EIIP_encoding=",EIIP_encoding.shape)
    return EIIP_encoding
    

def toMMI(inFasta,outCSV,A_type="DNA"):
    
    NA = 'ACGT'
    if A_type=="RNA":
        NA="ACGU"
    
    dinucleotide_list = [a1 + a2 for a1 in NA for a2 in NA]
    trinucleotide_list = [a1 + a2 + a3 for a1 in NA for a2 in NA for a3 in NA]
    dinucleotide_dict = {}
    trinucleotide_dict = {}
    
    for elem in dinucleotide_list:
        dinucleotide_dict[''.join(sorted(elem))] = 0
    for elem in trinucleotide_list:
        trinucleotide_dict[''.join(sorted(elem))] = 0

    encodings = []
    header = ['PID', 'label']
    header += ['MMI_%s' % elem for elem in sorted(dinucleotide_dict.keys())]
    header += ['MMI_%s' % elem for elem in sorted(trinucleotide_dict.keys())]
    
    fasta_list=read_fasta(inFasta)

    for i in fasta_list:
        name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
        code = [name, label]
        
        f1_dict = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0,
            
        }
        
        if A_type=="RNA":
        
            f1_dict = {
            'A': 0,
            'C': 0,
            'G': 0,
            
            'U': 0,
        }
        
        
        f2_dict = dinucleotide_dict.copy()
        f3_dict = trinucleotide_dict.copy()

        for elem in sequence:
            if elem in f1_dict:
                f1_dict[elem] += 1
        for key in f1_dict:
            f1_dict[key] /= len(sequence)

        for i in range(len(sequence) - 1):
            if ''.join(sorted(sequence[i: i + 2])) in f2_dict:
                f2_dict[''.join(sorted(sequence[i: i + 2]))] += 1
        for key in f2_dict:
            f2_dict[key] /= (len(sequence) - 1)

        for i in range(len(sequence) - 2):
            if ''.join(sorted(sequence[i: i + 3])) in f3_dict:
                f3_dict[''.join(sorted(sequence[i: i + 3]))] += 1
        for key in f3_dict:
            f3_dict[key] /= (len(sequence) - 2)

        for key in sorted(f2_dict.keys()):
            if f2_dict[key] != 0 and f1_dict[key[0]] * f1_dict[key[1]] != 0:
                code.append(f2_dict[key] * math.log(f2_dict[key] / (f1_dict[key[0]] * f1_dict[key[1]])))
            else:
                code.append(0)
                
        for key in sorted(f3_dict.keys()):
            element_1 = 0
            element_2 = 0
            element_3 = 0
            if f2_dict[key[0:2]] != 0 and f1_dict[key[0]] * f1_dict[key[1]] != 0:
                element_1 = f2_dict[key[0:2]] * math.log(f2_dict[key[0:2]] / (f1_dict[key[0]] * f1_dict[key[1]]))
            if f2_dict[key[0] + key[2]] != 0 and f1_dict[key[2]] != 0:
                element_2 = (f2_dict[key[0] + key[2]] / f1_dict[key[2]]) * math.log(
                    f2_dict[key[0] + key[2]] / f1_dict[key[2]])
            if f2_dict[key[1:3]] != 0 and f3_dict[key] / f2_dict[key[1:3]] != 0:
                element_3 = (f3_dict[key] / f2_dict[key[1:3]]) * math.log(f3_dict[key] / f2_dict[key[1:3]])
            code.append(element_1 + element_2 - element_3)
        encodings.append(code)
        
    encoding_array = np.array(encodings, dtype=str)
    MMI_encoding=pd.DataFrame(encoding_array)
    MMI_encoding.columns=header
    MMI_encoding.to_csv(outCSV,index=False)
    print("MMI is DONE. MMI_encoding=",MMI_encoding.shape)
    return MMI_encoding
    
    
    
def toNMBroto(inFasta,outCSV,nlag=2,A_type="DNA"):

    encoding_array = np.array([])
    file_name=""
    property_name=""
    
    if A_type == 'DNA':
        file_name = './util_data/didnaPhyche.data'
        diDNAPhyChe=pd.read_csv("./util_data/diDNAPhyche.txt",sep="\t",index_col=0,header=0)
        property_name = diDNAPhyChe.index
    else:
        file_name = './util_data/dirnaPhyche.data'
        diRNAPhyChe=pd.read_csv("./util_data/diRNAPhyche.txt",sep="\t",index_col=0,header=0)
        property_name = diRNAPhyChe.index
        
         
    try:
        data_file = os.path.split(os.path.realpath(__file__))[0] + file_name
        with open(data_file, 'rb') as handle:
            property_dict = pickle.load(handle)
            
    except Exception as e:
        error_msg = 'Could not find the physicochemical properties file.'
        return False
        
    nlag = nlag

    # value normalization
    for p_name in property_name:
        tmp = np.array(property_dict[p_name], dtype=float)
        pmean = np.average(tmp)
        pstd = np.std(tmp)
        property_dict[p_name] = [(elem - pmean) / pstd for elem in tmp]

    base = 'ACGT'
    if A_type=="RNA":
        base = 'ACGU'
        
    encodings = []
    header = ['PID', 'Label']
    for p_name in property_name:
        for d in range(1, nlag + 1):
            header.append(p_name + '.lag' + str(d))
     

    AADict = {}
    AA_list = [aa1 + aa2 for aa1 in base for aa2 in base]
    for i in range(len(AA_list)):
        AADict[AA_list[i]] = i
    
    fasta_list=read_fasta(inFasta)
    for elem in fasta_list:
        name, sequence, label = elem[0], re.sub('-', '', elem[1]), str(elem[2])
        code = [name, label]
        N = len(sequence) - 1
        for p_name in property_name:
            for d in range(1, nlag + 1):
                try:
                    if N > nlag:
                        atsd = sum([float(property_dict[p_name][AADict[sequence[j: j+2]]]) * float(property_dict[p_name][AADict[sequence[j+d: j+d+2]]]) for j in range(N-d)]) / (N - d)
                    else:
                        atsd = 0
                except Exception as e:
                    atsd = 0
                code.append(atsd)
        encodings.append(code)

    
    encoding_array = np.array(encodings, dtype=str)
    NMBroto_encoding=pd.DataFrame(encoding_array)
    NMBroto_encoding.columns=header
    NMBroto_encoding.to_csv(outCSV,index=False)
    print("NMBroto is DONE. NMBroto_encoding=",NMBroto_encoding.shape)
    return NMBroto_encoding

def toMoran(inFasta,outCSV,nlag=2,A_type="DNA"):

    encoding_array = np.array([])

    if A_type == 'DNA':
        file_name = './util_data/didnaPhyche.data'
        diDNAPhyChe=pd.read_csv("./util_data/diDNAPhyche.txt",sep="\t",index_col=0,header=0)
        property_name = diDNAPhyChe.index
    else:
        file_name = './util_data/dirnaPhyche.data'
        diRNAPhyChe=pd.read_csv("./util_data/diRNAPhyche.txt",sep="\t",index_col=0,header=0)
        property_name = diRNAPhyChe.index
        
                   
    with open(file_name, 'rb') as handle:
        property_dict = pickle.load(handle)
            
   
    nlag =  nlag 

    # value normalization
    for p_name in property_name:
         
        tmp = np.array(property_dict[p_name], dtype=float)
        pmean = np.average(tmp)
        pstd = np.std(tmp)
        property_dict[p_name] = [(elem - pmean) / pstd for elem in tmp]

    base = 'ACGT'
    if A_type=="RNA":
        base="ACGU"
        
    encodings = []
    header = ['PID', 'Label']
    for p_name in property_name:
        for d in range(1, nlag + 1):
            header.append(p_name + '.lag' + str(d))
 
   
    
    AADict = {}
    AA_list = [aa1 + aa2 for aa1 in base for aa2 in base]
    for i in range(len(AA_list)):
        AADict[AA_list[i]] = i
    
    fasta_list=read_fasta(inFasta)

    for elem in fasta_list:
        name, sequence, label = elem[0], re.sub('-', '', elem[1]), str(elem[2])
        code = [name, label]
        N = len(sequence) - 1
        for p_name in property_name:
            
            pmean = sum([property_dict[p_name][AADict[sequence[i: i+2]]] for i in range(N)]) / N
            for d in range(1, nlag + 1):
                try:
                    Idup = sum([(property_dict[p_name][AADict[sequence[j: j+2]]] - pmean) * (property_dict[p_name][AADict[sequence[j+d: j+d+2]]] - pmean) for j in range(N-d)]) / (N - d)
                    Iddown = sum([(property_dict[p_name][AADict[sequence[j: j+2]]] - pmean) ** 2 for j in range(N-d)]) / N
                    code.append(Idup/Iddown)
                except Exception as e:
                    code.append(0)
        encodings.append(code)

    encoding_array = np.array(encodings, dtype=str)
    Moran_encoding=pd.DataFrame(encoding_array)
    Moran_encoding.columns=header
    Moran_encoding.to_csv(outCSV,index=False)
    print("Moran is DONE. NMBroto_encoding=",Moran_encoding.shape)
    return Moran_encoding

def toGeary(inFasta,outCSV,nlag=2,A_type="DNA"):
    encoding_array = np.array([])

    if A_type == 'DNA':
        file_name = './util_data/didnaPhyche.data'
        diDNAPhyChe=pd.read_csv("./util_data/diDNAPhyche.txt",sep="\t",index_col=0,header=0)
        property_name = diDNAPhyChe.index
    else:
        file_name = './util_data/dirnaPhyche.data'
        diRNAPhyChe=pd.read_csv("./util_data/diRNAPhyche.txt",sep="\t",index_col=0,header=0)
        property_name = diRNAPhyChe.index
        
    
    
    with open(file_name, 'rb') as handle:
        property_dict = pickle.load(handle)
        
    nlag =  nlag 
    
    

    # value normalization
    for p_name in property_name:
        tmp = np.array(property_dict[p_name], dtype=float)
        pmean = np.average(tmp)
        pstd = np.std(tmp)
        property_dict[p_name] = [(elem - pmean) / pstd for elem in tmp]

    base = 'ACGT'
    if A_type=="RNA":
        base="ACGU"
    encodings = []
    header = ['PID', 'Label']
    for p_name in property_name:
        for d in range(1, nlag + 1):
            header.append(p_name + '.lag' + str(d))
    encodings.append(header)

    AADict = {}
    AA_list = [aa1 + aa2 for aa1 in base for aa2 in base]
    for i in range(len(AA_list)):
        AADict[AA_list[i]] = i
    
    
    fasta_list=read_fasta(inFasta)
    for elem in fasta_list:
        name, sequence, label = elem[0], re.sub('-', '', elem[1]),str(elem[2])
         
        code = [name, label]
        
        N = len(sequence) - 1
        for p_name in property_name:
            pmean = sum([property_dict[p_name][AADict[sequence[i: i + 2]]] for i in range(N)]) / N
            for d in range(1, nlag + 1):
                try:
                    Cdup = sum([(property_dict[p_name][AADict[sequence[j: j+2]]] - property_dict[p_name][AADict[sequence[j+d: j+d+2]]]) ** 2 for j in range(N - d)]) / (2 * (N - d))
                    Cddown = sum([(property_dict[p_name][AADict[sequence[j: j+2]]] - pmean) ** 2 for j in range(N - d)]) / (N - 1)
                    code.append(Cdup / Cddown)
                except Exception as e:
                    code.append(0)
        encodings.append(code)
        
    encoding_array = np.array(encodings, dtype=str)
    Geary_encoding=pd.DataFrame(encoding_array)
    Geary_encoding.columns=header
    Geary_encoding.to_csv(outCSV,index=False)
    print("Geary is DONE. Geary_encoding=",Geary_encoding.shape)
    return Geary_encoding
    


'''
自相关性
'''
def generatePropertyPairs( myPropertyName):
    pairs = []
    for i in range(len(myPropertyName)):
        for j in range(i + 1, len(myPropertyName)):
            pairs.append([myPropertyName[i], myPropertyName[j]])
            pairs.append([myPropertyName[j], myPropertyName[i]])
    return pairs

myDNADiIndex = {
            'AA': 0, 'AC': 1, 'AG': 2, 'AT': 3,
            'CA': 4, 'CC': 5, 'CG': 6, 'CT': 7,
            'GA': 8, 'GC': 9, 'GG': 10, 'GT': 11,
            'TA': 12, 'TC': 13, 'TG': 14, 'TT': 15
        }

myRNADiIndex = {
    'AA': 0, 'AC': 1, 'AG': 2, 'AU': 3,
    'CA': 4, 'CC': 5, 'CG': 6, 'CU': 7,
    'GA': 8, 'GC': 9, 'GG': 10, 'GU': 11,
    'UA': 12, 'UC': 13, 'UG': 14, 'UU': 15
}

didna_list = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content', 'A-philicity',
              'Propeller twist', 'Duplex stability:(freeenergy)',
              'Duplex tability(disruptenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
              'Stabilising energy of Z-DNA', 'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH',
              'Breslauer_dS', 'Electron_interaction', 'Hartman_trans_free_energy', 'Helix-Coil_transition',
              'Ivanov_BA_transition', 'Lisser_BZ_transition', 'Polar_interaction', 'SantaLucia_dG',
              'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Stacking_energy',
              'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS', 'Watson-Crick_interaction', 'Twist', 'Tilt', 'Roll',
              'Shift', 'Slide', 'Rise',
              'Clash Strength', 'Roll_roll', 'Twist stiffness', 'Tilt stiffness', 'Shift_rise',
              'Adenine content', 'Direction', 'Twist_shift', 'Enthalpy1', 'Twist_twist', 'Roll_shift',
              'Shift_slide', 'Shift2', 'Tilt3', 'Tilt1', 'Tilt4', 'Tilt2', 'Slide (DNA-protein complex)1',
              'Tilt_shift', 'Twist_tilt', 'Twist (DNA-protein complex)1', 'Tilt_rise', 'Roll_rise',
              'Stacking energy', 'Stacking energy1', 'Stacking energy2', 'Stacking energy3', 'Propeller Twist',
              'Roll11', 'Rise (DNA-protein complex)', 'Tilt_tilt', 'Roll4', 'Roll2', 'Roll3', 'Roll1',
              'Minor Groove Size', 'GC content', 'Slide_slide', 'Enthalpy', 'Shift_shift', 'Slide stiffness',
              'Melting Temperature1', 'Flexibility_slide', 'Minor Groove Distance',
              'Rise (DNA-protein complex)1', 'Tilt (DNA-protein complex)', 'Guanine content',
              'Roll (DNA-protein complex)1', 'Entropy', 'Cytosine content', 'Major Groove Size', 'Twist_rise',
              'Major Groove Distance', 'Twist (DNA-protein complex)', 'Purine (AG) content',
              'Melting Temperature', 'Free energy', 'Tilt_slide', 'Major Groove Width', 'Major Groove Depth',
              'Wedge', 'Free energy8', 'Free energy6', 'Free energy7', 'Free energy4', 'Free energy5',
              'Free energy2', 'Free energy3', 'Free energy1', 'Twist_roll', 'Shift (DNA-protein complex)',
              'Rise_rise', 'Flexibility_shift', 'Shift (DNA-protein complex)1', 'Thymine content', 'Slide_rise',
              'Tilt_roll', 'Tip', 'Keto (GT) content', 'Roll stiffness', 'Minor Groove Width', 'Inclination',
              'Entropy1', 'Roll_slide', 'Slide (DNA-protein complex)', 'Twist1', 'Twist3', 'Twist2', 'Twist5',
              'Twist4', 'Twist7', 'Twist6', 'Tilt (DNA-protein complex)1', 'Twist_slide', 'Minor Groove Depth',
              'Roll (DNA-protein complex)', 'Rise2', 'Persistance Length', 'Rise3', 'Shift stiffness',
              'Probability contacting nucleosome core', 'Mobility to bend towards major groove', 'Slide3',
              'Slide2', 'Slide1', 'Shift1', 'Bend', 'Rise1', 'Rise stiffness',
              'Mobility to bend towards minor groove']
              
dirna_list = ['Slide (RNA)', 'Adenine content', 'Hydrophilicity (RNA)', 'Tilt (RNA)', 'Stacking energy (RNA)',
              'Twist (RNA)', 'Entropy (RNA)', 'Roll (RNA)', 'Purine (AG) content', 'Hydrophilicity (RNA)1',
              'Enthalpy (RNA)1', 'GC content', 'Entropy (RNA)1', 'Rise (RNA)', 'Free energy (RNA)',
              'Keto (GT) content', 'Free energy (RNA)1', 'Enthalpy (RNA)', 'Guanine content', 'Shift (RNA)',
              'Cytosine content', 'Thymine content']

myDict = {
    'DAC': {'DNA': didna_list, 'RNA': dirna_list},
    'DCC': {'DNA': didna_list, 'RNA': dirna_list},
    'DACC': {'DNA': didna_list, 'RNA': dirna_list}
    }
    
myDataFile = {
    'DAC': {'DNA': './util_data/didnaPhyche.data', 'RNA': './util_data/dirnaPhyche.data'},
    'DCC': {'DNA': './util_data/didnaPhyche.data', 'RNA': './util_data/dirnaPhyche.data'},
    'DACC': {'DNA': './util_data/didnaPhyche.data', 'RNA': './util_data/dirnaPhyche.data'},
    }
    
def toDACvector(inFasta,outCSV, nlag=2 ,A_type="DNA"):
     
    if A_type=="DNA"  :
        myPropertyName=myDict["DAC"]["DNA"]
        myPropertyValue_file=myDataFile["DAC"]["DNA"]
        with open(myPropertyValue_file, 'rb') as handle:
            myPropertyValue = pickle.load(handle)
        myDiIndex = myDNADiIndex
        
    elif A_type=="RNA"  :
        myPropertyName=myDict["DAC"]["RNA"]
        myPropertyValue_file=myDataFile["DAC"]["RNA"]
        with open(myPropertyValue_file, 'rb') as handle:
            myPropertyValue = pickle.load(handle)
        print(myPropertyValue)
        myDiIndex = myRNADiIndex
     
    
    fastas = read_fasta(inFasta)
    
    lag =  nlag 
    
    kmer=2
    encodings = []
    myIndex = myDiIndex
    
    header = ['PID', 'Label']
    for p in myPropertyName:
        for l in range(1, lag + 1):
            header.append('%s.lag%d' % (p, l))
    
    

    for i in fastas:
        # print(i[0])
        name, sequence, label = i[0], re.sub('-', '', i[1]), str(i[2])
        code = [name, label]

        for p in myPropertyName:
            meanValue = 0
            # for j in range(len(sequence) - kmer):
            for j in range(len(sequence) - kmer + 1):
                  
                sub_seq=sequence[j: j + kmer]
                idx=myIndex[sub_seq]
                #print(type(p),p,type(sub_seq),sub_seq,type(idx),idx)
                #print(myPropertyValue[p][idx])
                
                #sss
                meanValue = meanValue + float(myPropertyValue[p][idx])
                #meanValue = meanValue + float(myPropertyValue[p][myIndex[sequence[j: j + kmer]]])
            # meanValue = meanValue / (len(sequence) - kmer)
            meanValue = meanValue / (len(sequence) - kmer + 1)

            for l in range(1, lag + 1):
                acValue = 0
                for j in range(len(sequence) - kmer - l + 1):
                    # acValue = acValue + (float(myPropertyValue[p][myIndex[sequence[j: j+kmer]]]) - meanValue) * (float(myPropertyValue[p][myIndex[sequence[j+l:j+l+kmer]]]))
                    acValue = acValue + (float(myPropertyValue[p][myIndex[sequence[j: j + kmer]]]) - meanValue) * (
                            float(myPropertyValue[p][myIndex[sequence[j + l:j + l + kmer]]]) - meanValue)
                acValue = acValue / (len(sequence) - kmer - l + 1)
                # print(acValue)
                code.append(acValue)
        encodings.append(code)
    
    encoding_array = np.array(encodings, dtype=str)
    DACvector_encoding=pd.DataFrame(encoding_array)
    DACvector_encoding.columns=header
    DACvector_encoding.to_csv(outCSV,index=False)
    print("DACvector is DONE. DACvector_encoding=",DACvector_encoding.shape)
    return DACvector_encoding
    
    
    
def toDCCvector(inFasta,outCSV, nlag=2 ,A_type="DNA"):
 
    kmer=2
    lag=nlag
    if A_type=="DNA"  :
        myPropertyName=myDict["DAC"]["DNA"]
        myPropertyValue_file=myDataFile["DAC"]["DNA"]
        with open(myPropertyValue_file, 'rb') as handle:
            myPropertyValue = pickle.load(handle)
        myDiIndex = myDNADiIndex
        
    elif A_type=="RNA"  :
        myPropertyName=myDict["DAC"]["RNA"]
        myPropertyValue_file=myDataFile["DAC"]["RNA"]
        with open(myPropertyValue_file, 'rb') as handle:
            myPropertyValue = pickle.load(handle)
        print(myPropertyValue)
        myDiIndex = myRNADiIndex
        
    encodings = []
    myIndex = myDiIndex if kmer == 2 else myTriIndex
    
    propertyPairs = generatePropertyPairs(myPropertyName)
    header = ['PID', 'Label'] + [n[0] + '-' + n[1] + '-lag.' + str(l) for n in propertyPairs for l in range(1, lag + 1)]
     
    fastas=read_fasta(inFasta)
    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), str(i[2])
        code = [name, label]

        for pair in propertyPairs:
            meanP1 = 0
            meanP2 = 0
            # for j in range(len(sequence) - kmer):
            for j in range(len(sequence) - kmer + 1):
                meanP1 = meanP1 + float(myPropertyValue[pair[0]][myIndex[sequence[j: j + kmer]]])
                meanP2 = meanP2 + float(myPropertyValue[pair[1]][myIndex[sequence[j: j + kmer]]])
            # meanP1 = meanP1 / (len(sequence) - kmer)
            # meanP2 = meanP2 / (len(sequence) - kmer)
            meanP1 = meanP1 / (len(sequence) - kmer + 1)
            meanP2 = meanP2 / (len(sequence) - kmer + 1)

            for l in range(1, lag + 1):
                ccValue = 0
                for j in range(len(sequence) - kmer - l + 1):
                    ccValue = ccValue + (
                            float(myPropertyValue[pair[0]][myIndex[sequence[j: j + kmer]]]) - meanP1) * (
                                    float(
                                        myPropertyValue[pair[1]][myIndex[sequence[j + l:j + l + kmer]]]) - meanP2)
                ccValue = ccValue / (len(sequence) - kmer - l + 1)
                code.append(ccValue)
        encodings.append(code)
     
    encoding_array = np.array(encodings, dtype=str)
    DCCvector_encoding=pd.DataFrame(encoding_array)
    DCCvector_encoding.columns=header
    DCCvector_encoding.to_csv(outCSV,index=False)
    print("DCCvector is DONE. DCCvector_encoding=",DCCvector_encoding.shape)
    return DCCvector_encoding
    
    

def toDACCvector(inFasta,outCSV, nlag=2 ,A_type="DNA"):

    kmer=2
    lag=nlag
    if A_type=="DNA"  :
        myPropertyName=myDict["DAC"]["DNA"]
        myPropertyValue_file=myDataFile["DAC"]["DNA"]
        with open(myPropertyValue_file, 'rb') as handle:
            myPropertyValue = pickle.load(handle)
        myDiIndex = myDNADiIndex
        
    elif A_type=="RNA"  :
        myPropertyName=myDict["DAC"]["RNA"]
        myPropertyValue_file=myDataFile["DAC"]["RNA"]
        with open(myPropertyValue_file, 'rb') as handle:
            myPropertyValue = pickle.load(handle)
        print(myPropertyValue)
        myDiIndex = myRNADiIndex
        
    encodings = []
    myIndex = myDiIndex if kmer == 2 else myTriIndex
    if len(myPropertyName) < 2:
        print('Error: two or more property are needed for cross covariance (i.e. DCC and TCC) descriptors')
        sys.exit(1)

    header = ['PID', 'Label']
    for p in myPropertyName:
        for l in range(1, lag + 1):
            header.append('%s.lag%d' %(p, l))
            
    propertyPairs = generatePropertyPairs(myPropertyName)
    header = header + [n[0] + '-' + n[1] + '-lag.' + str(l) for n in propertyPairs for l in range(1, lag + 1)]
    
    
    fastas=read_fasta(inFasta)
    for i in fastas:
        name, sequence, label = i[0], re.sub('-', '', i[1]), str(i[2])
        code = [name, label]
        ## Auto covariance
        for p in myPropertyName:
            meanValue = 0
            # for j in range(len(sequence) - kmer):
            for j in range(len(sequence) - kmer + 1):
                meanValue = meanValue + float(myPropertyValue[p][myIndex[sequence[j: j + kmer]]])
            # meanValue = meanValue / (len(sequence) - kmer)
            meanValue = meanValue / (len(sequence) - kmer + 1)

            for l in range(1, lag + 1):
                acValue = 0
                for j in range(len(sequence) - kmer - l + 1):
                    # acValue = acValue + (float(myPropertyValue[p][myIndex[sequence[j: j+kmer]]]) - meanValue) * (float(myPropertyValue[p][myIndex[sequence[j+l:j+l+kmer]]]))
                    acValue = acValue + (float(myPropertyValue[p][myIndex[sequence[j: j + kmer]]]) - meanValue) * (
                        float(myPropertyValue[p][myIndex[sequence[j + l:j + l + kmer]]]) - meanValue)
                acValue = acValue / (len(sequence) - kmer - l + 1)
                # print(acValue)
                code.append(acValue)

        ## Cross covariance
        for pair in propertyPairs:
            meanP1 = 0
            meanP2 = 0
            #for j in range(len(sequence) - kmer):
            for j in range(len(sequence) - kmer + 1):
                meanP1 = meanP1 + float(myPropertyValue[pair[0]][myIndex[sequence[j: j+kmer]]])
                meanP2 = meanP2 + float(myPropertyValue[pair[1]][myIndex[sequence[j: j+kmer]]])
            #meanP1 = meanP1 / (len(sequence) - kmer)
            #meanP2 = meanP2 / (len(sequence) - kmer)
            meanP1 = meanP1 / (len(sequence) - kmer + 1)
            meanP2 = meanP2 / (len(sequence) - kmer + 1)

            for l in range(1, lag + 1):
                ccValue = 0
                for j in range(len(sequence) - kmer - l + 1):
                    ccValue = ccValue + (float(myPropertyValue[pair[0]][myIndex[sequence[j: j + kmer]]]) - meanP1) * (
                        float(myPropertyValue[pair[1]][myIndex[sequence[j + l:j + l + kmer]]]) - meanP2)
                ccValue = ccValue / (len(sequence) - kmer - l + 1)
                code.append(ccValue)
        encodings.append(code)
    
    encoding_array = np.array(encodings, dtype=str)
    DACCvector_encoding=pd.DataFrame(encoding_array)
    DACCvector_encoding.columns=header
    DACCvector_encoding.to_csv(outCSV,index=False)
    print("DACCvector is DONE. DACCvector_encoding=",DACCvector_encoding.shape)
    return DACCvector_encoding
    
    
     
    