# Co-occurrence Structures Visualizer
# Author: Zhenyu Guo
# Requirements:
# Bio [pip install biopython]
# reportlab [pip install reportlab]
# rlPyCairo [pip install rlPyCairo]

import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
# from reportlab.lib import colors
# from reportlab.lib.units import cm

Program_Dir_Path = os.path.dirname(os.path.abspath(__file__))
Fasta_Fname = f"{Program_Dir_Path}/figS1/RiskSeqs.ERR3209754.contigs.fasta"
Outname = f"{Program_Dir_Path}/figS1/RiskSeqs.Visualization.ERR3209754.contigs.pdf"

Mark_Color_Dict = {'1': "red", '2': "green", '3': "blue"}

def readFastaFile(fasta_fname: str):
    return {record.description: str(record.seq) for record in SeqIO.parse(fasta_fname, 'fasta')}

def drawSeqs(fasta_dict: dict, outname: str, seqnum=10, seq_minlen=1000, seq_maxlen=float("inf"), gene_minnum=3, lable_size=8):
    gdd = GenomeDiagram.Diagram("Co-ocurrence Structure Visualization")
    seqname_list = list(fasta_dict.keys())
    
    current_seqnum = 0
    for seqname in seqname_list:
        contig_name_tmpls = seqname.split('| ')
        contig_len = len(fasta_dict[seqname])
        contig_name, gene_info_ls = " {}, {} bp".format(contig_name_tmpls[0], contig_len), contig_name_tmpls[1].split('; ')
        if contig_len >= seq_minlen and contig_len <= seq_maxlen and len(gene_info_ls) >= gene_minnum:
            gdt_features = gdd.new_track(1, name=contig_name, start=0, end=contig_len, greytrack=True, greytrack_labels=1, greytrack_fontsize=8)
            gds_features_set = gdt_features.new_set()
            
            for gene_info in gene_info_ls:
                gene_tmpls = gene_info.split(',')
                gene, start_raw, end_raw = gene_tmpls[0], int(gene_tmpls[1]), int(gene_tmpls[2])
                if start_raw > end_raw:
                    strand, start, end = -1, end_raw, start_raw
                else:
                    strand, start, end = +1, start_raw, end_raw
                
                gene_mark = gene[0]
                gene = gene[1:]
                color = Mark_Color_Dict[gene_mark]
                if gene[:2] in ('Tn', 'In'):
                    gene = gene.split('-')[0]
                
                feature = SeqFeature(FeatureLocation(start, end, strand=strand))
                gds_features_set.add_feature(feature, name=gene, label=True, label_size=lable_size, label_position="start", label_angle=30,
                                        sigil="ARROW", color=color, arrowshaft_height=0.4)
            
            current_seqnum += 1
        
        if current_seqnum == seqnum:
            break
    
    gdd.draw(format="linear", pagesize="A3", fragments=1)
    gdd.write(outname, "pdf")

if __name__ == "__main__":
    drawSeqs(readFastaFile(Fasta_Fname), Outname)
