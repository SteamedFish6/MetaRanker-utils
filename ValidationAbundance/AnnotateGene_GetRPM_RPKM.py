# Annotate Genes and Get Abundance (RPM & RPKM) of Genes
# Author: Zhenyu Guo

# Input: Non-redundant Gene Sequences and Depth of Genes
# Output: RPM & RPKM of Genes

# Main process:
# 1.annotate genes using blast;
# 2.remove redundant gene sequences;
# 3.count gene depth from reads.

# Required tools (need to be installed in system):
# blastn (ncbi-blast+)

import os
import time
import subprocess
import numpy as np
import pandas as pd
from Bio import SeqIO


# input and output path
Sample_Gene_Path = "AbundValidation/GenePredict/GeneFasta"
Sample_Depth_Path = "AbundValidation/GenePredict/GeneDepth"
Out_Path = "AbundValidation/GenePredict"

def collectFile(path: str, isfolder=False, ftype='', fullpath=False) -> list:
    file_list = os.listdir(path)
    collect_list = []
    for file_name in file_list:
        full_pathname = os.path.join(path, file_name)
        if isfolder == True:
            if os.path.isdir(full_pathname):
                if fullpath == False:
                    collect_list.append(file_name)
                else:
                    collect_list.append(full_pathname)
        else:
            if file_name[-(len(ftype)+1):] == '.{}'.format(ftype):
                if fullpath == False:
                    collect_list.append(file_name)
                else:
                    collect_list.append(full_pathname)
    collect_list = sorted(collect_list)
    return collect_list

Query_File_List = []
Query_Depthfile_List = []
for Samplename in collectFile(Sample_Gene_Path, isfolder=True):
    Sample_Gene_File = os.path.join(Sample_Gene_Path, Samplename, f"{Samplename}.CDS.cdhit.fna")
    Sample_Depth_File = os.path.join(Sample_Depth_Path, Samplename, f"{Samplename}.genes.depth.tsv")
    Query_File_List.append(Sample_Gene_File)
    Query_Depthfile_List.append(Sample_Depth_File)

# BLAST config
Seq_Minlen = 0
Blast_Method = "blastn" #"blastp"
NUM_THREAD = 32
BLAST_EVALUE = 0.0001
PERC_IDENTITY = 80 #not used in blastp
is_Cover_Old_File = False

# database config, when use_inner_refgene == False, use refgene_abund_list_dict_nucl & refgene_abund_list_dict_prot instead
Db_Path = "ranker_db/db_fasta"
Db_List = ["CARD", "MGE", "VFDB"]
Colinear_Min_NonCross_Len = 75

# filter criteria
Filter_Criteria = 'identity>=85 & cover_len>=75'

# backpasted source files, which corresponding to query_file_list
Seqname_Lendict_Fname = "ranker_db/ranker_blastdb_seqname_length.tsv"

Category_Fname_Dict = {
                       "CARD": "ranker_db/aro_index.tsv",
                       "VFDB": "ranker_db/VFs.tsv",
                       }
is_Category_Sepline = True

# abundance calculating config
is_Calc_Abund = True
is_Calc_RPM = True
is_Calc_RPKM = True

Refgene_Abund_List = [[] for _ in range(len(Query_Depthfile_List))]
Sample_List = [[] for _ in range(len(Query_Depthfile_List))]


def PrepareBlastDB(dbpath: str, dblist: str, blast_method: str, is_cover=False) -> None:
    if blast_method == "blastn":
        extend_name = "nhr"
        dbtype = "nucl"
    if blast_method == "blastp":
        extend_name = "phr"
        dbtype = "prot"
        dblist = ["{}_pro".format(db) for db in dblist]
    print("Preparing BLAST Database...")
    for dbname in dblist:
        dbfile = os.path.join(dbpath, "{}.fasta".format(dbname))
        checkCallCMD("makeblastdb -dbtype {} -parse_seqids -in {}".format(dbtype, dbfile), "{}.{}".format(dbfile, extend_name), force=is_cover)
    print("BLAST Database is Ready.")

class Query:
    def __init__(self, query_fname: str, blastmethod: str, dbpath: str, dblist: list, 
                 depth_df: pd.DataFrame=pd.DataFrame(), samplenames: list=None, refgene_abund: list=None):
        self.query_fname = query_fname
        self.blastmethod = blastmethod
        self.dbpath = dbpath
        self.dblist = dblist
        self.depth_df = depth_df
        if samplenames:
            self.depth_df.columns = samplenames
        self.samplenames = list(depth_df.columns)
        self.total_reads_num = self.depth_df.sum(axis=0)
        self.refgene_abund = refgene_abund
        self.query_tag = os.path.basename(query_fname).rsplit('.', 1)[0]
        
        self.M8_fdict = {} #key: db, value: M8_fname
        self.M8pp_fdict = {} #key: db, value: M8_fname
        self.rank_fdict = {} #key: samplename ('_' is contig with no depth), value: Risk_df_outname
        
        self.abs_abund_fdict = {} #key: db, value: abund_fname
    
    def FastaPreprocess(self, outpath: str, minlen=0) -> None:
        print("Checking input fasta file...")
        record_list = readFastaFile(self.query_fname, to_which='list')
        if minlen <= 0:
            ncontig = len(record_list)
            print("{} sequences found in {}".format(ncontig, self.query_fname))
        else:
            new_record_list = []
            for record in record_list:
                if len(record.seq) >= minlen:
                    new_record_list.append(record)
            outname = os.path.join(outpath, "filter{}_{}".format(minlen, os.path.basename(self.query_fname)))
            SeqIO.write(new_record_list, outname, "fasta")
            print("{} was filtered by >= {} bp/aa, new fasta is {}".format(self.query_fname, minlen, outname))
            self.query_fname = outname
            ncontig = len(new_record_list)
            print("{} sequences found in {}".format(ncontig, self.query_fname))
        
        if self.depth_df.empty == False:
            print("Checking input depth file...")
            for samplename in self.samplenames:
                sample_depth = self.depth_df[samplename]
                ncontig_sample = len(sample_depth[sample_depth > 0])
                print("{} sequences found in sample: {}".format(ncontig_sample, samplename))
    
    def Blast(self, outpath: str, is_cover=False) -> None:
        if self.blastmethod == 'blastn':
            perc_identity_arg = "-perc_identity {} ".format(PERC_IDENTITY)
            a_dblist = self.dblist.copy()
        elif self.blastmethod == 'blastp':
            perc_identity_arg = ""
            a_dblist = ["{}_pro".format(db) for db in self.dblist]
        for dbname, a_dbname in zip(self.dblist, a_dblist):
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), end=' ')
            db_fname = os.path.join(self.dbpath, "{}.fasta".format(a_dbname))
            print("Processing {}: {}, in db: {}".format(self.blastmethod, self.query_fname, db_fname))
            outname = os.path.join(outpath, "out_{}_{}_{}.tsv".format(self.blastmethod, self.query_tag, dbname))
            checkCallCMD("{} -task {} -query {} -db {} -out {} -outfmt \"6\" -evalue {} {}-num_threads {}".format(self.blastmethod, self.blastmethod, self.query_fname, db_fname, outname, BLAST_EVALUE, perc_identity_arg, NUM_THREAD), outname, force=is_cover)
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), end=' ')
            print("Finished {}: {}, in db: {}".format(self.blastmethod, self.query_fname, db_fname))
            self.M8_fdict[dbname] = outname
    
    def M8Preprocess(self, outpath: str, lendict: dict, is_pos_max=False, is_group_max=False, criteria: str=None, sortby='score', noncross_len=75) -> None:
        print("Handling BLAST result...")
        column_index = ['query_ID', 'target_ID', 'identity', 'cover_len', 'not_match', 'gap', 'start_query', 'end_query', 'start_target', 'end_target', 'evalue', 'score']
        for dbname in self.dblist:
            M8file = self.M8_fdict[dbname]
            
            # FilterM8
            sheet = pd.read_csv(M8file, sep='\t', header=None, names=column_index)
            outname = os.path.join(outpath, "agm_bp_sfa_"+os.path.basename(M8file))
            if criteria != None:
                sheet = sheet.query(criteria)
            sheet = sheet.sort_values(by=sortby, ascending=False)
            
            if is_group_max == True: # sort, ascending=False, retain the elements appear for the first time
                outname = os.path.join(outpath, "agm_bp_sfg_"+os.path.basename(M8file))
                newsheet = pd.DataFrame(columns=column_index)
                groupmax_name_dict = {}
                for i in sheet.index:
                    row = sheet.loc[i]
                    groupmax_name = row['query_ID']
                    if groupmax_name not in groupmax_name_dict:
                        groupmax_name_dict[groupmax_name] = 1
                        newsheet.loc[i] = row
                sheet = newsheet
            
            elif is_pos_max == True: # for same contig, remove duplicates by start and end positions
                outname = os.path.join(outpath, "agm_bp_sfp_"+os.path.basename(M8file))
                newsheet = pd.DataFrame(columns=column_index)
                all_contigs = sheet['query_ID'].unique()
                for contig_ID in all_contigs:
                    contig_df = sheet[sheet['query_ID']==contig_ID]
                    ##shape
                    if contig_df.shape[0] == 1:
                        newsheet.loc[contig_df.index[0]] = contig_df.iloc[0]
                    else:  # record covered range with (min, max)
                        contig_df = contig_df.sort_values(by=sortby, ascending=False)
                        newsheet.loc[contig_df.index[0]] = contig_df.iloc[0]
                        s1, e1 = sorted((contig_df.iloc[0]['start_query'], contig_df.iloc[0]['end_query']))
                        for i in contig_df.index[1:]:
                            contig_row = contig_df.loc[i]
                            s2, e2 = sorted((contig_row['start_query'], contig_row['end_query']))
                            if s1 - s2 >= noncross_len:
                                newsheet.loc[i] = contig_row
                                s1 = s2
                            elif e2 - e1 >= noncross_len:
                                newsheet.loc[i] = contig_row
                                e1 = e2
                sheet = newsheet
            
            # Name/LengthBackpaste & AddGeneName
            target_id_list = []
            target_len_list = []
            gene_name_list = []
            old_target_ids = sheet['target_ID']
            for masked_target_ID in old_target_ids:
                target_id_list.append(lendict[masked_target_ID][0]) #targetID_original
                target_len_list.append(eval(lendict[masked_target_ID][1])) #target_len
                gene_name_list.append(lendict[masked_target_ID][2])
            sheet['target_ID'] = target_id_list
            del target_id_list
            sheet.insert(2, 'gene_name', gene_name_list)
            del gene_name_list
            sheet.insert(3, 'target_len', target_len_list)
            del target_len_list
            
            # DepthBackpaste
            if self.depth_df.empty == False:
                M8_depth_df = pd.DataFrame(index=sheet.index, columns=self.samplenames)
                for i in sheet.index:
                    query_ID = sheet.loc[i, 'query_ID']
                    try:
                        depth = self.depth_df.loc[query_ID]
                    except:
                        pass
                    else:
                        M8_depth_df.loc[i] = depth
                M8_depth_df = M8_depth_df.astype(np.float64)
                M8_depth_df.fillna(0.0, inplace=True)
                for samplename in self.samplenames:
                    sheet[samplename] = M8_depth_df[samplename]
            
            sheet.to_csv(outname, index=False, sep='\t')
            self.M8pp_fdict[dbname] = outname
    
    def AddCateName(self, outpath: str, cate_dict: dict, sepline=False) -> None:
        print("Adding category annotation to BLAST result...")
        cate_list_dict = {"CARD": ['AMR Gene Family', 'Drug Class', 'Resistance Mechanism'],
                          "VFDB": ['VF_name', 'VF_bacteria', 'VF_category'],}
        checksep_dict_dict = {"CARD": {'AMR Gene Family': False, 'Drug Class': True, 'Resistance Mechanism': True},
                              "VFDB": {'VF_name': False, 'VF_bacteria': False, 'VF_category': False},}
        
        for dbname in self.dblist:
            if dbname in cate_list_dict:
                cate_names = cate_list_dict[dbname]
                cate_names_num = len(cate_names)
                M8_fname = self.M8pp_fdict[dbname]
                sheet = pd.read_csv(M8_fname, sep='\t')
                col_index = sheet.columns
                
                for cate_i, cate_type in enumerate(cate_names):
                    outname = os.path.join(outpath, f"{os.path.basename(M8_fname).rsplit('.', 1)[0]}_{cate_type}.tsv")
                    newsheet = pd.DataFrame(columns=col_index)
                    checksep = checksep_dict_dict[dbname][cate_type]
                    for i in range(len(sheet.index)):
                        row = sheet.iloc[i].copy()
                        cate_num = extractCateNum(row['target_ID'], dbname)
                        try:
                            cate_name = renameCateName(cate_dict[cate_num][cate_i], dbname, cate_type)
                        except:
                            cate_name = 'undefined'
                        
                        if sepline == True and checksep == True and ';' in cate_name:
                            ARnameLS = cate_name.split(';')
                        else:
                            ARnameLS = [cate_name]
                        for cate_name in ARnameLS:
                            row['gene_name'] = cate_name
                            newsheet = pd.concat([newsheet, row.to_frame().T], ignore_index=True)
                    
                    self.M8pp_fdict[f"{dbname}_{cate_type}"] = outname
                    newsheet.to_csv(outname, index=False, sep='\t')
                
                add_cols = [[] for idx in range(cate_names_num)]
                outname_all = os.path.join(outpath, f"{os.path.basename(M8_fname).rsplit('.', 1)[0]}_alltype.tsv")
                
                for i in range(len(sheet.index)):
                    row = sheet.iloc[i]
                    cate_num = extractCateNum(row['target_ID'], dbname)
                    for idx in range(cate_names_num):
                        try:
                            cate_name = renameCateName(cate_dict[cate_num][idx], dbname, cate_names[idx])
                        except:
                            cate_name = 'undefined'
                        add_cols[idx].append(cate_name)
                    
                for idx in range(cate_names_num):
                    sheet.insert(idx+3, cate_names[idx], add_cols[idx])
                
                sheet.to_csv(outname_all, index=False, sep='\t')
            
            elif dbname == "MGE":
                M8_fname = self.M8pp_fdict[dbname]
                outname = os.path.join(outpath, f"{os.path.basename(M8_fname).rsplit('.', 1)[0]}_MGE_category.tsv")
                sheet = pd.read_csv(M8_fname, sep='\t')
                col_index = sheet.columns
                gene_names = sheet['gene_name']
                cate_name_list = []
                for genename in gene_names:
                    cate_name_list.append(renameCateName(genename, "MGE"))
                sheet['gene_name'] = cate_name_list
                
                self.M8pp_fdict["MGE_category"] = outname
                sheet.to_csv(outname, index=False, sep='\t')
    
    def CalcAbsAbund(self, outpath: str, is_sorted=True) -> None:
        print("Calculating absolute abundance...")
        for dbname in self.M8pp_fdict:
            fname = self.M8pp_fdict[dbname]
            outname = os.path.join(outpath, f"AbsAbund.{self.query_tag}.{dbname}.tsv")
            sheet = pd.read_csv(fname, sep='\t')
            col_index = ['gene'] + self.samplenames
            abund_sheet = pd.DataFrame(columns=col_index)
            abund_sheet = abund_sheet.set_index('gene')
            
            coverage_col = sheet['cover_len'] / sheet['target_len'] # covered length / reference length
            reads_num_sheet = sheet.iloc[:, -len(self.samplenames):]
            reads_num_sheet.astype(np.float64)
            
            for i in range(len(sheet.index)):
                reads_num_series = reads_num_sheet.iloc[i]
                gene_name = sheet.loc[sheet.index[i], 'gene_name']
                coverage = coverage_col[i]
                abund_series = reads_num_series * coverage
                if gene_name not in abund_sheet.index:
                    abund_sheet.loc[gene_name] = abund_series
                else:
                    abund_sheet.loc[gene_name] = abund_sheet.loc[gene_name] + abund_series
            
            if is_sorted:
                abund_sheet['row_sum'] = abund_sheet.sum(axis=1)
                abund_sheet = abund_sheet.sort_values(by='row_sum', ascending=False)
            
            abund_sheet.to_csv(outname, columns=self.samplenames, sep='\t')
            self.abs_abund_fdict[dbname] = outname
    
    
    def CalcRPM_RPKM(self, outpath: str, calc_type='RPM', is_sorted=True) -> None:
        if self.depth_df.empty == True:
            pass
        else:
            print("Calculating {} abundance...".format(calc_type))
            total_reads_num = self.total_reads_num
            for dbname in self.M8pp_fdict:
                fname = self.M8pp_fdict[dbname]
                outname = os.path.join(outpath, f"{calc_type}.{self.query_tag}.{dbname}.tsv")
                sheet = pd.read_csv(fname, sep='\t')
                col_index = ['gene'] + self.samplenames
                abund_sheet = pd.DataFrame(columns=col_index)
                abund_sheet = abund_sheet.set_index('gene')
                
                reads_num_sheet = sheet.iloc[:, -len(self.samplenames):]
                reads_num_sheet.astype(np.float64)
                
                for i in range(len(sheet.index)):
                    reads_num_series = reads_num_sheet.iloc[i]
                    gene_name = sheet.loc[sheet.index[i], 'gene_name']
                    target_len = sheet.loc[sheet.index[i], 'target_len']
                    abund_series = reads_num_series * 1000000 / total_reads_num # 'RPM'
                    if calc_type == 'RPKM':
                        abund_series = abund_series * 1000 / target_len# 'RPKM'
                    
                    if gene_name not in abund_sheet.index:
                        abund_sheet.loc[gene_name] = abund_series
                    else:
                        abund_sheet.loc[gene_name] = abund_sheet.loc[gene_name] + abund_series
                
                if is_sorted:
                    abund_sheet['row_sum'] = abund_sheet.sum(axis=1)
                    abund_sheet = abund_sheet.sort_values(by='row_sum', ascending=False)
                abund_sheet.to_csv(outname, columns=self.samplenames, sep='\t')
    
    def CalcRelativeAbund(self, outpath: str) -> None:
        print("Calculating relative abundance...")
        for dbname in self.abs_abund_fdict:
            fname = self.abs_abund_fdict[dbname]
            outname = os.path.join(outpath, f"RelAbund.{self.query_tag}.{dbname}.tsv")
            sheet = pd.read_csv(fname, sep='\t', index_col='gene')
            sum_abund_list = []
            for samplename in self.samplenames: # sum by column, calculate total abundance
                sum_abund_list.append(sheet[samplename].sum())
            for sample_id, samplename in enumerate(self.samplenames):
                sheet[samplename] = sheet[samplename] / sum_abund_list[sample_id]
            sheet.to_csv(outname, sep='\t')


def checkCallCMD(cmd: str, outname: str, force=False) -> None:
    '''If outname not exists, then call cmd
    '''
    if os.path.exists(outname) == False or force == True:
        subprocess.check_call(cmd, shell=True)

def makeDir(dirparh: str, subname: str) -> str:
    make_path = os.path.join(dirparh, subname)
    if os.path.isdir(make_path) == False:
        os.mkdir(make_path)
    return make_path

def readSheetFile(fname: str, sep='\t', remove_header=True) -> list:
    tsv_sheet = []
    f1 = open(fname, 'r')
    for line in f1.readlines():
        tsv_line = line[:-1].split(sep)
        tsv_sheet.append(tsv_line)
    f1.close()
    if remove_header == True:
        return tsv_sheet[1:]# Remove titie
    else:
        return tsv_sheet

def loadLenDict(fname: str) -> dict:
    len_dict = {}
    len_sheet = readSheetFile(fname)
    for line in len_sheet:
        len_dict[line[0]] = [line[1], line[2], line[3]]
    return len_dict

def loadDepthDF(fname: str) -> pd.DataFrame:
    depth_df = pd.read_csv(fname, sep='\t', index_col=0)
    rowsum = np.sum(depth_df, axis=1)
    depth_df = depth_df.loc[rowsum > 0]
    return depth_df

def loadCateDict(fname_dict: dict) -> dict:
    cate_dict = {}
    for dbname in fname_dict:
        cate_sheet = readSheetFile(fname_dict[dbname])
        if dbname == "CARD":
            for line in cate_sheet:
                cate_dict[line[0]] = [line[8], line[9], line[10]]
        elif dbname == "VFDB":
            for line in cate_sheet:
                cate_dict[line[0]] = [line[1], line[3], line[5]]
    return cate_dict

def readFastaFile(fasta_fname: str, to_which='dict'):
    if to_which == 'dict':
        return {record.description: str(record.seq) for record in SeqIO.parse(fasta_fname, 'fasta')}
    elif to_which == 'list':
        return [record for record in SeqIO.parse(fasta_fname, "fasta")]

def extractCateNum(target_id: str, dbname: str) -> str:
    cate_num = target_id
    if dbname == "CARD":
        cate_num = "ARO:" + target_id.split(':')[1][:7]
    elif dbname == "VFDB":
        cate_num = 'VF' + target_id.split(' (VF')[1].split(') -')[0]
    elif dbname == "GREENGENES_16S":
        cate_num = target_id.split(' ', 1)[0]
    return cate_num

def renameCateName(cate_name: str, dbname: str, cate_type: str=None) -> str:
    if dbname == "CARD":
        if cate_type == "AMR Gene Family":
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


if __name__ == "__main__":
    
    PrepareBlastDB(Db_Path, Db_List, Blast_Method, is_cover=is_Cover_Old_File)
    Tmp_Path = makeDir(Out_Path, "temp")
    M8_Path = makeDir(Out_Path, "output_M8")
    M8pp_Path = makeDir(Out_Path, "preprocessed_M8")
    M8cate_Path = makeDir(M8pp_Path, "categorized_M8")
    
    Abund_Path = makeDir(Out_Path, "abundance")
    Abs_Abund_Path = makeDir(Abund_Path, "absolute")
    Relative_Abund_Path = makeDir(Abund_Path, "relative")
    RPM_Path = makeDir(Abund_Path, "RPM") if is_Calc_RPM else ''
    RPKM_Path = makeDir(Abund_Path, "RPKM") if is_Calc_RPKM else ''
    
    Len_Dict = loadLenDict(Seqname_Lendict_Fname)
    Cate_Dict = loadCateDict(Category_Fname_Dict)
    
    for i in range(len(Query_File_List)):
        Depth_Df = loadDepthDF(Query_Depthfile_List[i])
        Query_Obj = Query(Query_File_List[i], Blast_Method, Db_Path, Db_List, depth_df=Depth_Df, samplenames=Sample_List[i], refgene_abund=Refgene_Abund_List[i])
        Query_Obj.FastaPreprocess(minlen=Seq_Minlen, outpath=Tmp_Path)
        Query_Obj.Blast(outpath=M8_Path, is_cover=is_Cover_Old_File)
        Query_Obj.M8Preprocess(outpath=M8pp_Path, lendict=Len_Dict, is_group_max=False, is_pos_max=True, 
                                criteria=Filter_Criteria, noncross_len=Colinear_Min_NonCross_Len)
        Query_Obj.AddCateName(outpath=M8cate_Path, cate_dict=Cate_Dict, sepline=is_Category_Sepline)
        
        if is_Calc_Abund:
            Query_Obj.AddCateName(outpath=M8cate_Path, cate_dict=Cate_Dict, sepline=is_Category_Sepline)
            Query_Obj.CalcAbsAbund(outpath=Abs_Abund_Path, is_sorted=True)
            Query_Obj.CalcRelativeAbund(outpath=Relative_Abund_Path)
            if is_Calc_RPM:
                Query_Obj.CalcRPM_RPKM(RPM_Path, calc_type='RPM', is_sorted=True)
            if is_Calc_RPKM:
                Query_Obj.CalcRPM_RPKM(RPKM_Path, calc_type='RPKM', is_sorted=True)


