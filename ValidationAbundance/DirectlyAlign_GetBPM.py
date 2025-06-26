# Directly Align and Get Abundance (BPM)
# Author: Zhenyu Guo

# Input: Clean Reads
# Output: BPM of Genes

# Main process:
# 1.align reads to DB;
# 2.count gene depth from reads.

# Required tools (need to be installed in system):
# pigz
# bowtie2
# samtools

## Note: Please run `bowtie2 index` first!  

import os
import subprocess
import gzip
import pandas as pd


Sample_List_Fname = "AbundValidation/code/Sample_List.txt"
Reads_Path = "AbundValidation/ReadsFile"
Outpath = "AbundValidation/directly_align2"
Total_Collect_File = "AbundValidation/BPM_collection2.tsv"
Reads_Folder_List = []
with open(Sample_List_Fname) as sr:
    Samplename_List = [raw_Samplename.strip() for raw_Samplename in sr.readlines()]
for Samplename in Samplename_List:
    Sample_Reads_Folder = os.path.join(Reads_Path, Samplename)
    Reads_Folder_List.append(Sample_Reads_Folder)

NUM_THREAD = 32

card_name = "ranker_db/db_fasta/CARD.fasta"
mge_name = "ranker_db/db_fasta/MGE.fasta"
vfdb_name = "ranker_db/db_fasta/VFDB.fasta"
align_ls = [card_name, mge_name, vfdb_name]
seqname_length_fname = "ranker_db/ranker_blastdb_seqname_length.tsv"

def checkCallCMD(cmd: str, outname: str, is_cover_old=False) -> bool:
    '''If outname not exists, then call cmd
    '''
    if os.path.exists(outname) == False or is_cover_old == True:
        subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    else:
        return False

def makeDir(dirparh: str, subname: str) -> str:
    make_path = os.path.join(dirparh, subname)
    if os.path.isdir(make_path) == False:
        os.mkdir(make_path)
    return make_path

def countFastqReads_Bases(fname: str, is_gzipped=True):
    if is_gzipped:
        handle = gzip.open(fname, "rt")
    else:
        handle = open(fname, 'r')
    total_reads_num = 0
    double_seq_len = 0
    for line in handle:
        if line[0] == '@':
            total_reads_num += 1
        if line[0] not in ('@', '+', ' ', '\n'):
            double_seq_len += (len(line)-1)
    total_seq_len = double_seq_len // 2
    return total_reads_num, total_seq_len
def countFastqReads_Bases_pigzAWK(fname: str, is_gzipped=True):
    '''Try using pigz & awk to accelerate gzipped file reading and reads & bases counting'''
    awk_cmd = r"awk 'BEGIN {reads=0;bases=0} NR%4==1 && /^@/ {reads++} NR%4==2 {bases+=length} END {print reads,bases}'"
    if is_gzipped:
        cmd = "pigz -p {} -dck {} | {}".format(NUM_THREAD, fname, awk_cmd)
    else:
        cmd = "cat {} | {}".format(fname, awk_cmd)
    
    try:
        proc = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        print("Parsing (gzipped) reads file: failed to run shell command, trying default functions.")
        total_reads_num, total_seq_len = countFastqReads_Bases(fname, is_gzipped=is_gzipped)
    else:
        output = proc.stdout.strip()
        if not output:
            raise ValueError("Wrong format of fastq file")
        total_reads_num, total_seq_len = map(int, map(float, output.split())) #e.g., '5.43923e+09' must be converted to float first, then to int
    
    if total_reads_num == 0:
        raise ValueError("Num of total reads is 0, please check input reads file.")
    return total_reads_num, total_seq_len

if __name__ == "__main__":
    seqname_length_df = pd.read_csv(seqname_length_fname, sep='\t', index_col=0)
    seqname_length_dict = seqname_length_df['blastdbSeqLen'].to_dict()
    seqname_origin_dict = seqname_length_df['blastdbSeqID_original'].to_dict()
    
    for reads_folder in Reads_Folder_List:
        sample_name = os.path.basename(reads_folder)
        sample_path = makeDir(Outpath, sample_name)
        reads1 = os.path.join(reads_folder, "{}_1.clean.fastq.gz".format(sample_name))
        reads2 = os.path.join(reads_folder, "{}_2.clean.fastq.gz".format(sample_name))
        
        reads1_seqnum, reads1_basenum = countFastqReads_Bases_pigzAWK(reads1)
        reads2_seqnum, reads2_basenum = countFastqReads_Bases_pigzAWK(reads2)
        total_bases_num = reads1_basenum + reads2_basenum
        
        with open(Total_Collect_File, 'a') as tcf:
            tcf.write(f"{sample_name}")
        
        for align_target in align_ls:
            db_name = os.path.basename(align_target.split('.')[0])
            bam_fname = os.path.join(sample_path, "{}.{}.bam".format(sample_name, db_name))
            depth_fname = os.path.join(sample_path, "{}.{}.depth.tsv".format(sample_name, db_name))
            
            # cmd1 = "bwa mem -t {} -v 0 {} {} {} | samtools view -bS --threads {} | samtools sort --threads {} -o {}".format(NUM_THREAD, align_target, reads1, reads2, NUM_THREAD, NUM_THREAD, bam_fname)
            cmd1 = "bowtie2 --no-unal -p {} -x {} -1 {} -2 {} | samtools view -bS --threads {} | samtools sort --threads {} -o {}".format(NUM_THREAD, align_target.rsplit('.', 1)[0], reads1, reads2, NUM_THREAD, NUM_THREAD, bam_fname)
            checkCallCMD(cmd1, bam_fname)
            checkCallCMD("samtools index {}".format(bam_fname), depth_fname)
            cmd2 = "samtools depth -a {} | ".format(bam_fname) + r'''awk 'BEGIN {OFS="\t"} {sum[$1]+=$3; count[$1]++} END {for (gene in sum) print gene, sum[gene]/count[gene]}' ''' + "> {}".format(depth_fname)
            checkCallCMD(cmd2, depth_fname)
            
            if os.path.exists(bam_fname):
                os.remove(bam_fname)
                os.remove(bam_fname+'.bai')
        
            depth_df = pd.read_csv(depth_fname, sep='\t', header=None, index_col=0, names=['Depth'])
            depth_df.index.name = 'Element'
            bpm_dict = {}
            orig_dict = {}
            for idx in depth_df.index:
                depth = depth_df.loc[idx, 'Depth']
                bpm_dict[idx] = depth * seqname_length_dict[idx] / total_bases_num * 1e6
                orig_dict[idx] = seqname_origin_dict[idx]
            depth_df['BPM'] = bpm_dict
            depth_df.insert(1, 'OriginalName', orig_dict)
            
            bpm_sum = depth_df['BPM'].sum()
            with open(Total_Collect_File, 'a') as tcf:
                tcf.write(f"\t{bpm_sum}")
            
            bpm_fname = os.path.join(sample_path, "{}.{}.depth.bpm.tsv".format(sample_name, db_name))
            depth_df.to_csv(bpm_fname, sep='\t')
        
        with open(Total_Collect_File, 'a') as tcf:
            tcf.write(f"\n")
