import gzip
import random
from Bio import SeqIO

pos_r1 = "./SRR22519794/SRR22519794_1.clean.fastq.gz"
pos_r2 = "./SRR22519794/SRR22519794_2.clean.fastq.gz"
neg_r1 = "./insilico_reads/resource/QLFW/QLFW_1.clean.fq"
neg_r2 = "./insilico_reads/resource/QLFW/QLFW_2.clean.fq"
fasta_fname_list1 = ["./insilico_reads/resource/high_risk_seqs/ESKAPE.fasta",
                    "./insilico_reads/resource/high_risk_seqs/mcr.fasta",
                    "./insilico_reads/resource/high_risk_seqs/multi_drug.fasta",
]
fasta_fname_list2 = ["./insilico_reads/resource/low_risk_seqs/soil_bacteria_genome_complete.fa"]
seed_list = [1, 2, 3]
out_r_path = "./insilico_reads/resource/MockReads2"
rate_dict = {"25":0.25, "50":0.5, "75":0.75}

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

def readFastqFile(fname: str, is_gzipped=True, prefix='1'):
    if is_gzipped:
        handle = gzip.open(fname, "rt")
    else:
        handle = open(fname, 'r')
    
    fastq_dict = {}
    i = 0
    n = 0
    for line in handle:
        line = line.strip()
        i += 1
        if i % 4 == 1:
            read_name = line[1:]
        elif i % 4 == 2:
            seq = line
        elif i % 4 == 0:
            base_qual = line
            fastq_dict[f"{read_name}{prefix}"] = (seq, base_qual)
            n += 1
    
    return fastq_dict

def writeFastqFromDict(fastq_dict: dict, outname: str):
    # with gzip.open(outname, 'wt') as w:
    with open(outname, 'w') as w:
        for readname in fastq_dict:
            batch = f"@{readname}\n{fastq_dict[readname][0]}\n+\n{fastq_dict[readname][1]}\n"
            w.write(batch)

def readFastaFile(fasta_fname: str, to_which='dict'):
    if to_which == 'dict':
        return {record.description: str(record.seq) for record in SeqIO.parse(fasta_fname, 'fasta')}
    elif to_which == 'list':
        return [record for record in SeqIO.parse(fasta_fname, "fasta")]

def generateReads(pos_dict: dict, neg_dict: dict, pos_rate=1.0) -> dict:
    pos_num = int(len(pos_dict) * pos_rate)
    neg_num = int(len(neg_dict) * (1 - pos_rate))
    pos_keys = random.sample(list(pos_dict.keys()), pos_num)
    neg_keys = random.sample(list(neg_dict.keys()), neg_num)
    
    res_dict = {}
    for k in pos_keys:
        res_dict[k] = pos_dict[k]
    for k in neg_keys:
        res_dict[k] = neg_dict[k]
    
    return res_dict

def getReadsFromFasta(fasta_dict: dict, seqnum=6e7, read_len=150, read_vib=3, contig_len=5000, contig_vib=2000, contig_gap=50, min_contig_len=1000) -> dict:
    cut_dict = {}
    total_template_bases = 0
    
    for seqname in fasta_dict:
        seq = fasta_dict[seqname]
        start = 0
        while start < len(seq):
            segment_len = max(round(random.normalvariate(contig_len, contig_vib)), min_contig_len)
            segment_seq = seq[start: start+segment_len]
            start = start + segment_len + contig_gap
            segment_name = f"{seqname}#{start}:{start + segment_len}"
            if len(segment_seq) >= min_contig_len:
                cut_dict[segment_name] = segment_seq
                total_template_bases += segment_len
    
    reads_dict = {}
    for segment_name in cut_dict:
        segment_seq = cut_dict[segment_name]
        segment_len = len(segment_seq)
        sample_times = round(segment_len / total_template_bases * seqnum)
        for n in range(sample_times):
            cut_read_len = sorted([read_len-2*read_vib, round(random.normalvariate(read_len, read_vib)), read_len+2*read_vib])[1]
            max_start_pos = segment_len - cut_read_len
            start_pos = random.randrange(0, max_start_pos)
            read_seq = segment_seq[start_pos: start_pos+cut_read_len]
            read_name = f"{segment_name}#{n}"
            reads_dict[read_name] = (read_seq, "F"*cut_read_len)
    
    return reads_dict


if __name__ == "__main__":
    
    Pos1 = readFastqFile(pos_r1, prefix='1')
    Pos2 = readFastqFile(pos_r2, prefix='2')
    Pos1.update(Pos2)
    del Pos2
    print(f"Pos1:{len(Pos1)}")
    
    
    # Neg1 = readFastqFile(neg_r1, is_gzipped=False, prefix='1')
    # Neg2 = readFastqFile(neg_r2, is_gzipped=False, prefix='2')
    # Neg1.update(Neg2)
    # del Neg2
    # print(f"Neg1:{len(Neg1)}")
    
    for seed_num in seed_list:
        random.seed(seed_num)
        
        # Pos1_fasta_dict = {}
        # for fasta_fname in fasta_fname_list1:
        #     Pos_fasta_dict = readFastaFile(fasta_fname)
        #     Pos1_fasta_dict.update(Pos_fasta_dict)
        # print(f"Fasta:{len(Pos1_fasta_dict)}")
        # Pos1 = getReadsFromFasta(Pos1_fasta_dict, seqnum=6e7, read_len=150, read_vib=3, contig_len=5000, contig_vib=2000, contig_gap=50, min_contig_len=1000)
        # print(f"Pos1:{len(Pos1)}")
        
        Neg1_fasta_dict = {}
        for fasta_fname in fasta_fname_list2:
            Neg_fasta_dict = readFastaFile(fasta_fname)
            Neg1_fasta_dict.update(Neg_fasta_dict)
        print(f"Fasta:{len(Neg1_fasta_dict)}")
        Neg1 = getReadsFromFasta(Neg1_fasta_dict, seqnum=6e7, read_len=150, read_vib=3, contig_len=5000, contig_vib=2000, contig_gap=50, min_contig_len=1000)
        print(f"Neg1:{len(Neg1)}")
        
        # out_mark = f"HRG_LRG_{seed_num}"
        out_mark = f"SRR22519794_LRG_{seed_num}"
        # Outname = f"{out_r_path}/{out_mark}_{len(Pos1)}_100.fq"
        # writeFastqFromDict(Pos1, Outname)
        
        # Outname = f"{out_r_path}/{out_mark}_{len(Neg1)}_00.fq"
        # writeFastqFromDict(Neg1, Outname)
        
        for rate_name in rate_dict:
            rate = rate_dict[rate_name]
            Outdict = generateReads(Pos1, Neg1, rate)
            Outlen = len(Outdict)
            Outname = f"{out_r_path}/{out_mark}_{Outlen}_{rate_name}.fq"
            writeFastqFromDict(Outdict, Outname)
            del Outdict

        # del Pos1
        del Neg1
