# Summary High Risk ARGs, MGEs and VFs by counting BPM abundance and co-ocurrence frequence
# Author: Zhenyu Guo



import os
import pandas as pd

Blastdb_fname = "ranker_db/ranker_blastdb_seqname_length.tsv"
ARO_index_fname = "ranker_db/aro_index.tsv"
VFs_index_fname = "ranker_db/VFs.tsv"

RE_BPM_path = "metaranker_output/BPM"
RiskMatrix_path = "metaranker_output/risk_result/RiskMatrix"

Sample_List_fname = "SampleNames.txt" # Need to prepare first!
Out_fname = "HighRiskGeneSummary.tsv"

def buildSumDF(dbfname: str) -> pd.DataFrame:
    '''gene    dbtype    BPM    coocur
    '''
    db_df = pd.read_csv(dbfname, sep='\t', index_col=0)
    db_df.drop_duplicates('genename', inplace=True)
    tag_dict = {"CARD": 'ARG', "MGE": 'MGE', "VFDB": 'VF'}
    gene_list = db_df['genename'].to_list()
    dbtype_list = [tag_dict[tag.split('_', 1)[0]] for tag in db_df.index]
    
    db_cate_list1 = []
    db_cate_list2 = []
    db_cate_list3 = []
    for i in range(len(gene_list)):
        db_cate_list1.append(getCateNames(Cate_Dict, gene_list[i], dbtype_list[i], cate_no=0))
        db_cate_list2.append(getCateNames(Cate_Dict, gene_list[i], dbtype_list[i], cate_no=1))
        db_cate_list3.append(getCateNames(Cate_Dict, gene_list[i], dbtype_list[i], cate_no=2))
    
    new_df = pd.DataFrame({'dbtype': dbtype_list, 'category1': db_cate_list1, 
                           'category2': db_cate_list2, 'category3': db_cate_list3, 
                           'BPM': 0.0, 'coocur': 0.0, 'insample': 0.0}, index=gene_list)
    new_df.index.name = 'gene'
    
    return new_df

def loadSampleList(fname: str) -> list:
    spname_list = []
    with open(fname, 'r') as f:
        for line in f:
            spname_list.append(line.strip())
    return spname_list

def loadBPMDepth(fname: str) -> pd.DataFrame:
    '''RE    BPM    contig    gene
#     '''
    df = pd.read_csv(fname, sep='\t', index_col=0)
    REcompo_list = [RE_name.split('#') for RE_name in df.index]
    
    new_df = pd.DataFrame(index=df.index)
    new_df.index.name = 'RE'
    new_df['BPM'] = df['BPM']
    new_df['contig'] = [REcompo[0] for REcompo in REcompo_list]
    new_df['gene'] = [REcompo[1] for REcompo in REcompo_list]
    return new_df

def getCoocurFreq(fname: str) -> dict:
    '''#contig    ARG    MGE    VF
    contig  freq(with other REs)
    '''
    df = pd.read_csv(fname, sep=',', index_col=0)
    # df.columns = ['ARG', 'MGE', 'VF']
    freq_df = df.sum(axis=1) - 1
    # freq_df = df.sum(axis=1)
    freq_dict = freq_df.to_dict()
    return freq_dict


def readSheetFile(fname: str, sep='\t', remove_header=True) -> list:
    tsv_sheet = []
    f1 = open(fname, 'r')
    for line in f1.readlines():
        tsv_line = line[:-1].split(sep)
        tsv_sheet.append(tsv_line)
    f1.close()
    if remove_header == True:
        return tsv_sheet[1:] #Remove titie
    else:
        return tsv_sheet

def loadCateDict(fname_dict: dict) -> dict: #pandas is slower here
    cate_dict = {}
    for dbname in fname_dict:
        cate_sheet = readSheetFile(fname_dict[dbname])
        if dbname == "CARD":
            for line in cate_sheet:
                cate_dict[line[-1]] = [line[8], line[9], line[10]]
        elif dbname == "VFDB":
            for line in cate_sheet:
                cate_dict[line[1]] = [line[1], line[3], line[5]]
    return cate_dict

def renameCateName(cate_name: str, dbname: str, cate_no: int=None) -> str:
    if dbname == "ARG":
        if cate_no == 0:
            if 'antibiotic efflux pump' in cate_name:
                sp = cate_name.split('(')[1].split(')')[0]
                cate_name = sp + ' type drug efflux'
            elif 'beta-lactamase' in cate_name:
                cate_name = 'beta-lactam'
            elif 'tetracycline' in cate_name:
                cate_name = 'Tetracycline'
            elif len(cate_name) < 10:
                cate_name = 'Aminoglycoside'
            else:
                cate_name = 'Others'
        elif cate_no == 1:
            cate_name = cate_name.split(';')[0]
    
    elif dbname == "MGE":
        prefix = cate_name.split('_', 1)[0]
        # if prefix in ('plasmid', 'vir', 'proph'):
        #     cate_name = prefix #ACLAME
        if prefix[:2] == 'IS':
            cate_name = 'IS'
        elif prefix[:3] in ('Inc', 'rep', 'Col'):
            cate_name = 'plasmid'
        elif prefix[:4] == 'MITE':
            cate_name = 'IS'
        elif prefix[:2] in ('Tn', 'In', 'kappa'):
            cate_name = 'Tn'
        elif prefix[:9] == 'INTEGRALL':
            cate_name = 'INTEGRALL'
        elif prefix[0] == 'p':
            cate_name = 'plasmid'
        # else:
        #     cate_name = 'ICEberg'
    
    return cate_name

def getCateNames(cate_dict: dict, gene_name, db_name, cate_no: int=None):
    try:
        if db_name == "ARG":
            cate_name_raw = cate_dict[gene_name][cate_no]
            cate_name = renameCateName(cate_name_raw, db_name, cate_no)
        elif db_name == "MGE":
            cate_name = renameCateName(gene_name, db_name)
        elif db_name == "VF":
            cate_name = cate_dict[gene_name][cate_no]
        else:
            cate_name = "Unknown"
    except:
        cate_name = "Error"
    return cate_name
        

Sample_List = loadSampleList(Sample_List_fname)
Cate_Dict = loadCateDict({"CARD": ARO_index_fname, "VFDB": VFs_index_fname})

Sum_df = buildSumDF(Blastdb_fname)
for Sample_name in Sample_List:
    ARG_depth_fname = os.path.join(RE_BPM_path, f"BPM.{Sample_name}.CARD.tsv")
    MGE_depth_fname = os.path.join(RE_BPM_path, f"BPM.{Sample_name}.MGE.tsv")
    VF_depth_fname = os.path.join(RE_BPM_path, f"BPM.{Sample_name}.VFDB.tsv")
    RiskMatrix_fname = os.path.join(RiskMatrix_path, f"RiskMatrix_{Sample_name}.csv")
    
    if os.path.exists(ARG_depth_fname) and os.path.exists(MGE_depth_fname) and os.path.exists(VF_depth_fname):
        ARG_depth_df = loadBPMDepth(ARG_depth_fname)
        MGE_depth_df = loadBPMDepth(MGE_depth_fname)
        VF_depth_df = loadBPMDepth(VF_depth_fname)
        Sample_depth_df = pd.concat([ARG_depth_df, MGE_depth_df, VF_depth_df])
        Sample_depth_df['tmpRE'] = Sample_depth_df.index
        Sample_depth_df.drop_duplicates('tmpRE', inplace=True)
        Sample_depth_df.drop('tmpRE', axis=1, inplace=True)
        
        Sample_coocur_dict = getCoocurFreq(RiskMatrix_fname)
        
        for Risk_element in Sample_depth_df.index:
            Gene_name = Sample_depth_df.loc[Risk_element, 'gene']
            BPM_num = Sample_depth_df.loc[Risk_element, 'BPM']
            Contig_name = Sample_depth_df.loc[Risk_element, 'contig']
            Coocur_freq = Sample_coocur_dict[Contig_name]
            
            Sum_df.loc[Gene_name, 'BPM'] += BPM_num
            Sum_df.loc[Gene_name, 'coocur'] += Coocur_freq
        
        Gene_List = Sample_depth_df['gene'].unique().tolist()
        for Gene_name in Gene_List:
            Sum_df.loc[Gene_name, 'insample'] += 1
        
Sum_df.to_csv(Out_fname, sep='\t')

