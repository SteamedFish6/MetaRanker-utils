# DatabaseBuilder, build MetaRanker database from raw DB
# Author: Zhenyu Guo

# Main process:
# 1.remove redundancy in ARG, MGE, VF;
# 2.remove overlap ARG in VF;
# 3.remove overlap VF, ARG in MGE;
# 4.rename fasta sequences and extract gene name to `ranker_blastdb_seqname_length.tsv`

# Required tools (need to be installed in system):
# cd-hit
# blastn (ncbi-blast+)

import os
import subprocess
import pandas as pd

NUM_THREAD = 32
f_criteria_nucl = 'identity>=85 & cover_len>=75'
Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))

raw_ARG_fasta = f"{Program_Dir_Path}/raw/CARD.fasta"
raw_MGE_fasta = f"{Program_Dir_Path}/raw/MGE.fasta"
raw_VF_fasta = f"{Program_Dir_Path}/raw/VFDB.fasta"

cdhit_ARG_fasta = f"{Program_Dir_Path}/temp/fasta/cdhit.CARD.fasta"
cdhit_MGE_fasta = f"{Program_Dir_Path}/temp/fasta/cdhit.MGE.fasta"
cdhit_VF_fasta = f"{Program_Dir_Path}/temp/fasta/cdhit.VFDB.fasta"

renamed_cdhit_ARG_fasta = f"{Program_Dir_Path}/temp/fasta/renamed_cdhit.CARD.fasta" #!rename
renamed_cdhit_MGE_fasta = f"{Program_Dir_Path}/temp/fasta/renamed_cdhit.MGE.fasta"
renamed_cdhit_VF_fasta = f"{Program_Dir_Path}/temp/fasta/renamed_cdhit.VFDB.fasta"

blast_query_ARG = renamed_cdhit_ARG_fasta
blast_query_MGE = renamed_cdhit_MGE_fasta

blast_db_ARG = renamed_cdhit_ARG_fasta
blast_db_VF = renamed_cdhit_VF_fasta

blast_output_ARG_VF = f"{Program_Dir_Path}/temp/output_M8/out.blastn.cdhit_CARD.cdhit_VFDB.tsv"
blast_output_MGE_ARG = f"{Program_Dir_Path}/temp/output_M8/out.blastn.cdhit_MGE.cdhit_CARD.tsv"
blast_output_MGE_VF = f"{Program_Dir_Path}/temp/output_M8/out.blastn.cdhit_MGE.cdhit_VFDB.tsv"

filtered_blast_output_ARG_VF = f"{Program_Dir_Path}/temp/output_M8/sf_out.blastn.cdhit_CARD.cdhit_VFDB.tsv"
filtered_blast_output_MGE_ARG = f"{Program_Dir_Path}/temp/output_M8/sf_out.blastn.cdhit_MGE.cdhit_CARD.tsv"
filtered_blast_output_MGE_VF = f"{Program_Dir_Path}/temp/output_M8/sf_out.blastn.cdhit_MGE.cdhit_VFDB.tsv"

marked_ARG_fasta = f"{Program_Dir_Path}/temp/fasta/marked_cdhit_CARD.fasta"
marked_MGE_fasta = f"{Program_Dir_Path}/temp/fasta/marked_cdhit_MGE.fasta"

corrected_ARG_fasta = f"{Program_Dir_Path}/temp/corrected/corrected_cdhit_CARD.fasta"
corrected_MGE_fasta = f"{Program_Dir_Path}/temp/corrected/corrected_cdhit_MGE.fasta"
# corrected_VF_fasta = f"{Program_Dir_Path}/temp/corrected/corrected_cdhit_VFDB.fasta"
corrected_VF_fasta = renamed_cdhit_VF_fasta

final_ARG_fasta = f"{Program_Dir_Path}/final/CARD.fasta"
final_MGE_fasta = f"{Program_Dir_Path}/final/MGE.fasta"
final_VF_fasta = f"{Program_Dir_Path}/final/VFDB.fasta"
seqname_lendict_fname = f"{Program_Dir_Path}/final/ranker_blastdb_seqname_length.tsv"


def RemoveSpace(fname, outname):
    fasta_dict = readFastaFile(fname)
    new_fasta_dict = {}
    for seqname in fasta_dict:
        new_seqname = seqname.rsplit(' ', 1)[0]
        new_seqname = new_seqname.replace(' ', '_')
        new_seq = fasta_dict[seqname].replace(' ', '')
        new_fasta_dict[new_seqname] = new_seq
    writeFastaFromDict(new_fasta_dict, outname)

def FilterM8(M8file: str, outname: str, criteria: str=None) -> None:
    column_index = ['query_ID', 'target_ID', 'identity', 'cover_len', 'not_match', 'gap', 'start_query', 'end_query', 'start_target', 'end_target', 'evalue', 'score']
    # outname = os.path.join(filtered_path, "sfa_"+os.path.basename(M8file))
    sheet = pd.read_csv(M8file, sep='\t', header=None, names=column_index)
    if criteria != None:
        sheet = sheet.query(criteria)
    sheet = sheet.sort_values(by='score', ascending=False)
    
    # newsheet = pd.DataFrame(columns=column_index)
    # groupmax_name_list = []
    # for i in range(len(sheet.index)):
    #     groupmax_name = sheet.loc[sheet.index[i], 'query_ID']
    #     if groupmax_name not in groupmax_name_list:
    #         groupmax_name_list.append(groupmax_name)
    #         row = sheet.iloc[i]
    #         newsheet = pd.concat([newsheet, row.to_frame().T], ignore_index=True)
    # sheet = newsheet
    
    sheet.to_csv(outname, columns=column_index, index=False, header=True, mode='w', sep='\t')

def readFastaFile(fname: str) -> dict:
    f1 = open(fname, 'r')
    seq_dict = {}
    # num_seq = 0
    for line in f1:
        if line[0] == '>':
            seq_name = line[1:-1]
            seq_dict[seq_name] = ""
            # num_seq += 1
        elif line == '\n' or line[0] == '#':
            pass
        else:
            seq_dict[seq_name] = seq_dict[seq_name] + line[:-1]
    f1.close()
    #print(num_seq)
    return seq_dict

def writeFastaFromDict(fasta_dict: dict, outname: str) -> None:
    f2 = open(outname, 'w')
    for seqname in fasta_dict:
        f2.write(">{}\n".format(seqname))
        seq = fasta_dict[seqname].upper()
        
        sep_seqLS = []
        row_seqlen_max = 80
        while len(seq) >= row_seqlen_max:
            sep_seqLS.append(seq[:row_seqlen_max])
            seq = seq[row_seqlen_max:]
        if len(seq) > 0:
            sep_seqLS.append(seq)
        for sep_seq in sep_seqLS:
            f2.write(sep_seq+'\n')
    f2.close()

def MarkFastaSeqByPos(fasta_fname: str, M8_fname: str, outname: str) -> None:
    # column_index = ['query_ID', 'target_ID', 'identity', 'cover_len', 'not_match', 'gap', 'start_query', 'end_query', 'start_target', 'end_target', 'evalue', 'score']
    sheet = pd.read_csv(M8_fname, sep='\t')
    # sheet.set_index('query_ID', inplace=True)
    cut_pos_dict = {seqname:[] for seqname in sheet['query_ID']}
    for idx in sheet.index:
        query_ID = sheet.iloc[idx]['query_ID']
        start_pos = sheet.iloc[idx]['start_query']
        end_pos = sheet.iloc[idx]['end_query']
        cut_pos_dict[query_ID].append(sorted([int(start_pos), int(end_pos)]))
    
    fasta_dict = readFastaFile(fasta_fname)
    out_fasta_dict = {}
    
    changed_seqs = 0
    for seqname in fasta_dict:
        seq = fasta_dict[seqname]
        if seqname in cut_pos_dict:
            for cut_pos_list in cut_pos_dict[seqname]:
                start_pos = cut_pos_list[0] - 1
                end_pos = cut_pos_list[1]
                marked_seq = '-' * (end_pos - start_pos)
                seq = "".join([seq[:start_pos], marked_seq, seq[end_pos:]])
            changed_seqs += 1
        # if len(seq) >= minlen:
        out_fasta_dict[seqname] = seq
    writeFastaFromDict(out_fasta_dict, outname)

def CutMarkedSeq(fasta_fname: str, outname: str, minlen=75):
    fasta_dict = readFastaFile(fasta_fname)
    out_fasta_dict = {}
    for seqname in fasta_dict:
        seq = fasta_dict[seqname]
        seq = seq.replace('-', '')
        if len(seq) >= minlen:
            out_fasta_dict[seqname] = seq
    
    writeFastaFromDict(out_fasta_dict, outname)

def seprGeneName_new(old_name: str) -> str:
    if old_name[:3] == 'gb|': #ARG name
        gene_name = old_name.split(':')[1][8:].split('_')[0]
    elif old_name[:3] == 'VFG': #VFDB
        gene_name = old_name.split('_(')[1].split(')_')[0]
    elif old_name[:2] == 'IS':
        if old_name.count('_') == 2: #IS family
            # gene_name = 'IS_' + old_name.split('_')[2]
            gene_name = old_name.split('_')[0] #IS No
        else: #Tn accession number
            gene_name = 'Tn_' + old_name.split('-')[0].split('.')[0]
    elif old_name[:4] == 'gnl|': #INTEGRALL Integrase gene
        gene_name = 'INTEGRALL_' + old_name.split('|')[3]
    elif old_name[:3]  in ('Inc', 'rep', 'Col'):
        gene_name_compo = old_name.split('_')
        gene_name = "{}_{}".format(gene_name_compo[0], gene_name_compo[1])
    elif old_name[0] == 'p':
        gene_name_compo = old_name.split('_')
        gene_name = "{}_{}".format(gene_name_compo[0], gene_name_compo[1])
    # elif old_name[:5] == 'gene:': #aclame nucl
    #     gene_name1 = old_name.split(':')[1]
    #     if '# LocusTag: ' in old_name:
    #         gene_name2 = old_name.split('# LocusTag: ')[1].split(' ')[0]
    #     else:
    #         gene_name2 = 'others'
    #     gene_name = '{}_{}'.format(gene_name1, gene_name2)
    # elif old_name[:8] == 'protein:': #aclame prot
    #     gene_name1 = old_name.split(':')[1]
    #     if '# MgeName: ' in old_name:
    #         gene_name2 = old_name.split('# MgeName: ')[1].split(' ')[0]
    #     else:
    #         gene_name2 = 'others'
    #     gene_name = '{}_{}'.format(gene_name1, gene_name2)
    # elif old_name[:8] == 'ICEberg|': #ICE
    #     gene_name = 'ICE_' + old_name.split('|')[2]
    else:
        gene_name = old_name
    return gene_name

def createLenDictwithFasta(flist: list, outlist: list, mark_list: list, sheet_outname: str) -> None:
    # lensheet = []
    titlelist = ['blastdbSeqID_replaced', 'blastdbSeqID_original', 'blastdbSeqLen', 'genename']
    sheet = pd.DataFrame(columns=titlelist)
    nrow = 0
    
    for nf in range(len(flist)):
        filename = flist[nf]
        seq_dict = readFastaFile(filename)
        dbmark = mark_list[nf]
        outname_fasta = outlist[nf]
        
        fasta_dict = {}
        i = 0
        for seqname in seq_dict:
            seq = seq_dict[seqname]
            new_seqname = "{}_{}".format(dbmark, i)
            genename = seprGeneName_new(seqname)
            i += 1
            nrow += 1
            # lensheet.append([new_seqname, seqname, str(len(seq_dict[seqname])), genename])
            sheet.loc[nrow] = [new_seqname, seqname, len(seq), genename]
            fasta_dict[new_seqname] = seq
        writeFastaFromDict(fasta_dict, outname_fasta)
    
    # writeTSV(lensheet, sheet_outname, titlelist=titlelist)
    sheet.to_csv(sheet_outname, sep='\t', index=False)


if __name__ == "__main__":
    subprocess.call("cd-hit -i {} -o {} -c 0.85 -aS 0.85 -n 5 -d 0 -g 1 -M 0 -T {}".format(raw_ARG_fasta, cdhit_ARG_fasta, NUM_THREAD), shell=True)
    subprocess.call("cd-hit -i {} -o {} -c 0.85 -aS 0.85 -n 5 -d 0 -g 1 -M 0 -T {}".format(raw_MGE_fasta, cdhit_MGE_fasta, NUM_THREAD), shell=True)
    subprocess.call("cd-hit -i {} -o {} -c 0.85 -aS 0.85 -n 5 -d 0 -g 1 -M 0 -T {}".format(raw_VF_fasta, cdhit_VF_fasta, NUM_THREAD), shell=True)
    
    RemoveSpace(cdhit_ARG_fasta, renamed_cdhit_ARG_fasta)
    RemoveSpace(cdhit_MGE_fasta, renamed_cdhit_MGE_fasta)
    RemoveSpace(cdhit_VF_fasta, renamed_cdhit_VF_fasta)
    
    subprocess.call("makeblastdb -dbtype nucl -in {}".format(blast_db_ARG), shell=True)
    subprocess.call("makeblastdb -dbtype nucl -in {}".format(blast_db_VF), shell=True)
    
    subprocess.call("blastn -task blastn -query {} -db {} -out {} -outfmt \"6\" -evalue 0.0001 -perc_identity 80 -num_threads {}".format(blast_query_ARG, blast_db_VF, blast_output_ARG_VF, NUM_THREAD), shell=True)
    subprocess.call("blastn -task blastn -query {} -db {} -out {} -outfmt \"6\" -evalue 0.0001 -perc_identity 80 -num_threads {}".format(blast_query_MGE, blast_db_ARG, blast_output_MGE_ARG, NUM_THREAD), shell=True)
    subprocess.call("blastn -task blastn -query {} -db {} -out {} -outfmt \"6\" -evalue 0.0001 -perc_identity 80 -num_threads {}".format(blast_query_MGE, blast_db_VF, blast_output_MGE_VF, NUM_THREAD), shell=True)
    
    FilterM8(blast_output_ARG_VF, filtered_blast_output_ARG_VF, criteria=f_criteria_nucl)
    FilterM8(blast_output_MGE_ARG, filtered_blast_output_MGE_ARG, criteria=f_criteria_nucl)
    FilterM8(blast_output_MGE_VF, filtered_blast_output_MGE_VF, criteria=f_criteria_nucl)
    MarkFastaSeqByPos(blast_query_ARG, filtered_blast_output_ARG_VF, marked_ARG_fasta)
    MarkFastaSeqByPos(blast_query_MGE, filtered_blast_output_MGE_ARG, marked_MGE_fasta)
    MarkFastaSeqByPos(marked_MGE_fasta, filtered_blast_output_MGE_VF, marked_MGE_fasta)
    CutMarkedSeq(marked_ARG_fasta, corrected_ARG_fasta, minlen=75)
    CutMarkedSeq(marked_MGE_fasta, corrected_MGE_fasta, minlen=75)
    
    createLenDictwithFasta([corrected_ARG_fasta, corrected_MGE_fasta, corrected_VF_fasta], [final_ARG_fasta, final_MGE_fasta, final_VF_fasta], ["CARD", "MGE", "VFDB"], seqname_lendict_fname)
    
