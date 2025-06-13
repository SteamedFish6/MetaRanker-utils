# Predict Genes and Get Depth of Genes
# Author: Zhenyu Guo

# Input: Clean Reads and Assembled Contigs of Samples
# Output: Nucleotide and AminoAcid Gene Sequences and their Depth

# Main process:
# 1.predict genes from contigs;
# 2.remove redundant gene sequences;
# 3.count gene depth from reads.

# Required tools (need to be installed in system):
# prodigal
# cd-hit
# bwa
# samtools


import os
import subprocess
from multiprocessing import Pool
from Bio import SeqIO


THREADS = 32

Project_path = "/home/gzy/SRAdata/ValidateBatch3_20250609" #20250609

Reads_folder = "ReadsFile"
Assembly_folder = "Assembly" #Not full
GenePredict_folder = "GenePredict" #Not full
GeneFasta_folder = "GeneFasta" #Not full
GeneDepth_folder = "GeneDepth" #Not full
Temp_folder = "TempFile" #Not full


class SingleSample:
    
    def __init__(self, accession_ID: str, reads_path: str, contig_fname: str) -> None:
        self.accession_ID = accession_ID
        
        reads_fname1 = os.path.join(reads_path, "{}_1.clean.fastq.gz".format(self.accession_ID))
        reads_fname2 = os.path.join(reads_path, "{}_2.clean.fastq.gz".format(self.accession_ID))
        if os.path.exists(reads_fname1) and os.path.exists(reads_fname2):
            self.reads_type = 'pe'
            self.reads_fname1 = reads_fname1
            self.reads_fname2 = reads_fname2
        else:
            self.reads_type = 'se'
            self.reads_fname1 = os.path.join(reads_path, "{}.clean.fastq.gz".format(self.accession_ID))
            self.reads_fname2 = None
        
        self.contig_fname = contig_fname
        self.gene_nucl_fname = None
        self.gene_prot_fname = None
        self.cdhit_nucl_fname = None
        self.cdhit_prot_fname = None
        self.bam_fname = None
        self.depth_fname = None
    
    def prodigal(self, output_dir: str, temp_dir: str=None, threads=24):
        if temp_dir == None:
            temp_dir = output_dir
        
        runMultiProdigal(self.contig_fname, self.accession_ID, output_dir, temp_dir, threads=threads)
        ## rename seqs here
        gene_nucl_fname = os.path.join(output_dir, "{}.CDS.fna".format(self.accession_ID))
        gene_prot_fname = os.path.join(output_dir, "{}.prot.faa".format(self.accession_ID))
        renameFastaSeqs(gene_nucl_fname, gene_nucl_fname, self.accession_ID)
        renameFastaSeqs(gene_prot_fname, gene_prot_fname, self.accession_ID)
        self.gene_nucl_fname = gene_nucl_fname
        self.gene_prot_fname = gene_prot_fname
    
    def cdhit(self, output_dir: str, input_gene_nucl: str=None, input_gene_prot : str=None, threads=24):
        if input_gene_nucl != None and input_gene_prot != None:
            self.gene_nucl_fname = input_gene_nucl
            self.gene_prot_fname = input_gene_prot
        out_fname_nucl = os.path.join(output_dir, "{}.CDS.cdhit.fna".format(self.accession_ID))
        out_fname_prot = os.path.join(output_dir, "{}.prot.cdhit.faa".format(self.accession_ID))
        ## cd-hit the protein gene seqs, extract seqs name, and fetch nucleotide gene seqs
        subprocess.check_call("cd-hit -i {} -o {} -c 0.95 -aS 0.9 -n 5 -d 0 -g 1 -M 0 -T {}".format(self.gene_prot_fname, out_fname_prot, threads), shell=True)
        subprocess.check_call("rm {}.clstr".format(out_fname_prot), shell=True)
        self.cdhit_prot_fname = out_fname_prot
        
        ## fetch nucleotide gene seqs
        ## CDS and prot seqs' name are the same in prodigal output, and cd-hit will not rename (can extract from dict)
        pickSeqsFromDict(out_fname_prot, self.gene_nucl_fname, out_fname_nucl)
        self.cdhit_nucl_fname = out_fname_nucl
    
    def bwa(self, output_dir: str, threads=24):
        '''Align contigs.fa genes.fna or to reads.fq
        '''
        subprocess.check_call("bwa index {}".format(self.cdhit_nucl_fname), shell=True)
        out_bam_fname = os.path.join(output_dir, "{}.sorted.bam".format(self.accession_ID))
        if self.reads_type == 'pe':
            subprocess.check_call("bwa mem -t {} {} {} {} -v 2 | samtools view -bS --threads {} | samtools sort --threads {} -o {}".format(threads, self.cdhit_nucl_fname, self.reads_fname1, self.reads_fname2, threads, threads, out_bam_fname), shell=True)
        elif self.reads_type == 'se':
            subprocess.check_call("bwa mem -t {} {} {} -v 2 | samtools view -bS --threads {} | samtools sort --threads {} -o {}".format(threads, self.cdhit_nucl_fname, self.reads_fname1, threads, out_bam_fname), shell=True)
            
        subprocess.check_call("rm {}.amb".format(self.cdhit_nucl_fname), shell=True)
        subprocess.check_call("rm {}.ann".format(self.cdhit_nucl_fname), shell=True)
        subprocess.check_call("rm {}.bwt".format(self.cdhit_nucl_fname), shell=True)
        subprocess.check_call("rm {}.pac".format(self.cdhit_nucl_fname), shell=True)
        subprocess.check_call("rm {}.sa".format(self.cdhit_nucl_fname), shell=True)
        self.bam_fname = out_bam_fname
    
    def extract_depth2(self, output_dir: str):
        subprocess.check_call("samtools index {}".format(self.bam_fname), shell=True)
        sheet_outname = os.path.join(output_dir, "{}.genes.depth.tsv".format(self.accession_ID))
        subprocess.check_call("samtools depth -a {} | ".format(self.bam_fname) + '''awk 'BEGIN {OFS="\t"} {sum[$1]+=$3; count[$1]++} END {for (gene in sum) print gene, sum[gene]/count[gene]}' | ''' + "sed '1i Gene\tDepth' > {}".format(sheet_outname), shell=True)
        
        subprocess.check_call("rm {}".format(self.bam_fname), shell=True)
        subprocess.check_call("rm {}.bai".format(self.bam_fname), shell=True)
        self.depth_fname = sheet_outname


## FilePath Management
def makeDir(dirparh: str, subname: str) -> str:
    make_path = os.path.join(dirparh, subname)
    if os.path.isdir(make_path) == False:
        os.mkdir(make_path)
    return make_path

def collectFile(path: str, isfolder=False, ftype='tar.gz', fullpath=False) -> list:
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

def runProdigal(contig_fname: str, sample_name: str, out_dir: str):
    """Run prodigal.
    """
    out_pathname = os.path.join(out_dir, sample_name)
    subprocess.check_call("prodigal -i {} -o {}.gff -d {}.CDS.fna -a {}.prot.faa -f gff -p meta -m -q".format(contig_fname, out_pathname, out_pathname, out_pathname), shell=True)

def runMultiProdigal(contig_fname: str, sample_name: str, out_dir: str, temp_dir: str, threads=1):
    """Run prodigal in multiple threads.
    """
    out_pathname = os.path.join(out_dir, sample_name)
    if threads == 1:
        runProdigal(contig_fname, sample_name, out_dir)
    else:
        
        seprList2File(contig_fname, threads, temp_dir)
        pool = Pool(processes=threads)
        for i in range(1, threads+1):
            sep_contig_fname = os.path.join(temp_dir, f"{os.path.basename(contig_fname)}_{i}.fasta")
            sep_sample_name = f"{sample_name}_{i}"
            pool.apply_async(runProdigal, args=(sep_contig_fname, sep_sample_name, temp_dir))
        pool.close()
        pool.join()
        
        # delete the temporary files
        for i in range(1, threads+1):
            removeProdigalTemp(os.path.join(temp_dir, f"{os.path.basename(contig_fname)}_{i}.fasta"))
        
        # merge the files
        temp_pathname = os.path.join(temp_dir, sample_name)
        
        out_fname = f"{out_pathname}.CDS.fna"
        fname_list = [f"{temp_pathname}_{i}.CDS.fna" for i in range(1, threads+1)]
        catFiles(fname_list, out_fname)
        for fname in fname_list:
            removeProdigalTemp(fname)
        
        out_fname = f"{out_pathname}.prot.faa"
        fname_list = [f"{temp_pathname}_{i}.prot.faa" for i in range(1, threads+1)]
        catFiles(fname_list, out_fname)
        for fname in fname_list:
            removeProdigalTemp(fname)
        
        # out_fname = f"{out_pathname}.gff"
        fname_list = [f"{temp_pathname}_{i}.gff" for i in range(1, threads+1)]
        # catFiles(fname_list, out_fname)
        for fname in fname_list:
            removeProdigalTemp(fname)

def catFiles(flist: list, out_fname: str):
    with open(out_fname, 'w') as f_out:
        for in_fname in flist:
            if os.path.exists(in_fname):
                with open(in_fname, 'r') as f_in:
                    for line in f_in:
                        f_out.write(line)

def pickSeqsFromDict(prot_cdhit_fname: str, nucl_fname: str, out_nucl_fname: str):
    prot_cdhit_fasta_dict = readFastaFile(prot_cdhit_fname, to_which='dict')
    nucl_fasta_dict = readFastaFile(nucl_fname, to_which='dict')
    nucl_cdhit_fasta_dict = {}
    for seqname in nucl_fasta_dict:
        if seqname in prot_cdhit_fasta_dict:
            nucl_cdhit_fasta_dict[seqname] = nucl_fasta_dict[seqname]
    writeFastaFromDict(nucl_cdhit_fasta_dict, out_nucl_fname)

def removeProdigalTemp(fname: str) -> None:
    if os.path.exists(fname):
        os.remove(fname)

def renameFastaSeqs(raw_fasta_fname: str, new_fasta_fname: str, name_prefix: str, seq_min_len=0) -> None:
    raw_fasta_dict = readFastaFile(raw_fasta_fname, to_which='dict')
    new_fasta_dict = {}
    i = 0
    if seq_min_len == 0:
        for seqname in raw_fasta_dict:
            seq = raw_fasta_dict[seqname]
            i += 1
            new_fasta_dict["{}_{}".format(name_prefix, i)] = seq
    else:
        for seqname in raw_fasta_dict:
            seq = raw_fasta_dict[seqname]
            if len(seq) >= seq_min_len:
                i += 1
                new_fasta_dict["{}_{}".format(name_prefix, i)] = seq
    writeFastaFromDict(new_fasta_dict, new_fasta_fname)

def seprList2File(fasta_fname: str, num_split: int, temp_dir: str):
    """Split original fasta file into several fasta files.
    """
    record_list = readFastaFile(fasta_fname, to_which='list')
    sepr_len, resnum = divmod(len(record_list), num_split)
    if resnum != 0:
        sepr_len += 1
    
    # out_list = []
    num_files = 1
    for i in range(0, len(record_list), sepr_len):
        outname = os.path.join(temp_dir, f"{os.path.basename(fasta_fname)}_{num_files}.fasta")
        out_seq_list = record_list[i:i+sepr_len]
        SeqIO.write(out_seq_list, outname, "fasta")
        # out_list.append(outname)
        num_files += 1
    # return out_list


def readFastaFile(fasta_fname: str, to_which='dict'):
    if to_which == 'dict':
        return {record.description: str(record.seq) for record in SeqIO.parse(fasta_fname, 'fasta')}
    elif to_which == 'list':
        return [record for record in SeqIO.parse(fasta_fname, "fasta")]

def writeFastaFromDict(fasta_dict: dict, outname: str) -> None:
    f2 = open(outname, 'w')
    for item in fasta_dict:
        f2.write(">{}\n".format(item))
        seq = (fasta_dict[item]).upper()
        
        # if seq[-1] == '*':
        #     seq = seq[:-1]
        
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


if __name__ == "__main__":
    
    Reads_path = makeDir(Project_path, Reads_folder)
    Assembly_path = makeDir(Project_path, Assembly_folder)
    GenePredict_path = makeDir(Project_path, GenePredict_folder)
    GeneFasta_path = makeDir(GenePredict_path, GeneFasta_folder)
    GeneDepth_path = makeDir(GenePredict_path, GeneDepth_folder)
    Temp_path = makeDir(Project_path, Temp_folder)
    
    Samples_list = []
    for Samplename in collectFile(Reads_path, isfolder=True):
        Sample_Reads_path = os.path.join(Reads_path, Samplename)
        Sample_Contig_file = os.path.join(Assembly_path, Samplename, f"{Samplename}.contigs.fa")
        Sample_obj = SingleSample(Samplename, reads_path=Sample_Reads_path, contig_fname=Sample_Contig_file)
        Samples_list.append(Sample_obj)
    
    for Sample_obj in Samples_list:
        Sample_gene_fasta_path = makeDir(GeneFasta_path, Sample_obj.accession_ID)
        Sample_obj.prodigal(Sample_gene_fasta_path, temp_dir=Temp_path, threads=THREADS)
    
    for Sample_obj in Samples_list:
        Sample_gene_fasta_path = makeDir(GeneFasta_path, Sample_obj.accession_ID)
        Sample_obj.cdhit(Sample_gene_fasta_path, threads=THREADS)
        Sample_gene_depth_path = makeDir(GeneDepth_path, Sample_obj.accession_ID)
        Sample_obj.bwa(Sample_gene_depth_path, threads=THREADS)
        Sample_obj.extract_depth2(Sample_gene_depth_path)

