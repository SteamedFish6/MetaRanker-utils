import gzip
import random
# from Bio import SeqIO

# pos_r1 = "./SRR22519794/SRR22519794_1.clean.fastq.gz"
# pos_r2 = "./SRR22519794/SRR22519794_2.clean.fastq.gz"
# out_prefix = "SRR22519794sub"
# pos_r1 = "./QLFW/QLFW_1.clean.fastq.gz"
# pos_r2 = "./QLFW/QLFW_2.clean.fastq.gz"
# out_prefix = "QLFWsub"
# pos_r1 = "./ERR13925542/ERR13925542_1.clean.fastq.gz"
# pos_r2 = "./ERR13925542/ERR13925542_2.clean.fastq.gz"
# out_prefix = "ERR13925542sub"
# # pos_r1 = "./SRR17642938/SRR17642938_1.clean.fastq.gz"
# # pos_r2 = "./SRR17642938/SRR17642938_2.clean.fastq.gz"
# # out_prefix = "SRR17642938sub"
# pos_r1 = "./SRR18491219/SRR18491219_1.clean.fastq.gz"
# pos_r2 = "./SRR18491219/SRR18491219_2.clean.fastq.gz"
# out_prefix = "SRR18491219sub"
# pos_r1 = "./SRR26803261/SRR26803261_1.clean.fastq.gz"
# pos_r2 = "./SRR26803261/SRR26803261_2.clean.fastq.gz"
# out_prefix = "SRR26803261sub"
# #pos_r1 = "./SRR8208343/SRR8208343_1.clean.fastq.gz"
# #pos_r2 = "./SRR8208343/SRR8208343_2.clean.fastq.gz"
# #out_prefix = "SRR8208343sub"

# pos_r1 = "./SRR31677053/SRR31677053_1.clean.fastq.gz"
# pos_r2 = "./SRR31677053/SRR31677053_2.clean.fastq.gz"
# out_prefix = "SRR31677053sub"
# # pos_r1 = "./SRR18494403/SRR18494403_1.clean.fastq.gz"
# # pos_r2 = "./SRR18494403/SRR18494403_2.clean.fastq.gz"
# # out_prefix = "SRR18494403sub"
# pos_r1 = "./SRR23403438/SRR23403438_1.clean.fastq.gz"
# pos_r2 = "./SRR23403438/SRR23403438_2.clean.fastq.gz"
# out_prefix = "SRR23403438sub"
# pos_r1 = "./ERR4837146/ERR4837146_1.clean.fastq.gz"
# pos_r2 = "./ERR4837146/ERR4837146_2.clean.fastq.gz"
# out_prefix = "ERR4837146sub"
# pos_r1 = "./SRR23998354/SRR23998354_1.clean.fastq.gz"
# pos_r2 = "./SRR23998354/SRR23998354_2.clean.fastq.gz"
# out_prefix = "SRR23998354sub"
pos_r1 = "./SRR7686823/SRR7686823_1.clean.fastq.gz"
pos_r2 = "./SRR7686823/SRR7686823_2.clean.fastq.gz"
out_prefix = "SRR7686823sub"


is_gzipped_fq = True
seed_list = [1, 2, 3]
rate_dict = {'0005':0.005, '001':0.01, '002':0.02, '005':0.05, '010':0.1, '020':0.2, '040':0.4, '060':0.6, '080':0.8}
out_r_path = "./MockReadsSub2"

def readFastqFile(fname: str, is_gzipped=True): #, prefix='1'
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
            # fastq_dict[f"{read_name}{prefix}"] = (seq, base_qual)
            fastq_dict[read_name] = (seq, base_qual)
            n += 1
    
    return fastq_dict

def writeFastqFromDict(fastq_dict: dict, outname: str):
    # with gzip.open(outname, 'wt') as w:
    with open(outname, 'w') as w:
        for readname in fastq_dict:
            batch = f"@{readname}\n{fastq_dict[readname][0]}\n+\n{fastq_dict[readname][1]}\n"
            w.write(batch)

def pickReads2(fastq_dict1: dict, fastq_dict2: dict, rate=0.1):
    seq_num = int(len(fastq_dict1) * rate)
    seqno_all_list = [i for i in range(seq_num)]
    seqno_list = random.sample(seqno_all_list, seq_num)
    seqno_list.sort()
    
    fastq_dict_k1 = list(fastq_dict1.keys())
    seq_keys1 = [fastq_dict_k1[k] for k in seqno_list]
    res_dict1 = {k: fastq_dict1[k] for k in seq_keys1}
    fastq_dictk2 = list(fastq_dict2.keys())
    seq_keys2 = [fastq_dictk2[k] for k in seqno_list]
    res_dict2 = {k: fastq_dict2[k] for k in seq_keys2}
    
    return res_dict1, res_dict2

if __name__ == "__main__":
    Pos1 = readFastqFile(pos_r1, is_gzipped=is_gzipped_fq)
    print(f"Pos1:{len(Pos1)}")
    Pos2 = readFastqFile(pos_r2, is_gzipped=is_gzipped_fq)
    print(f"Pos2:{len(Pos2)}")
    for seed_num in seed_list:
        random.seed(seed_num)
        
        for rate_name in rate_dict:
            rate = rate_dict[rate_name]
            print(f"{seed_num}: {rate}")
            Outdict1, Outdict2 = pickReads2(Pos1, Pos2, rate)
            Outlen = len(Outdict1)
            Outname1 = f"{out_r_path}/{out_prefix}_{rate_name}_{seed_num}_{Outlen}_1.fq"
            Outname2 = f"{out_r_path}/{out_prefix}_{rate_name}_{seed_num}_{Outlen}_2.fq"
            writeFastqFromDict(Outdict1, Outname1)
            writeFastqFromDict(Outdict2, Outname2)
    
