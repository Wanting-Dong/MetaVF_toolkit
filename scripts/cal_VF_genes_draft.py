import sys
import os
class VFrecord(object):
    def __init__(self, line):
        fields = line.rstrip().split('\t')
        self.VF = fields[0]
        self.VFid = fields[1]
        self.category = fields[3]
        self.species = fields[4]
        self.species_all = fields[5]
        self.length = int(fields[8])
        self.ICE = fields[6]
        self.prophage = fields[7]
        self.plasmid = fields[9]

def create_VFid_dic(VFid_info_file):
    info_dict = {}
    with open(VFid_info_file,"r") as fileobject:
        next(fileobject)
        for line in fileobject:
            info_dict[VFrecord(line).VFid]=[VFrecord(line).species, VFrecord(line).plasmid, VFrecord(line).ICE, VFrecord(line).prophage,VFrecord(line).category]
    return info_dict


def cal_VFid_draft(VFid_info_file,blast_filer_file,output_file1,output_file2):
    info_dict = create_VFid_dic(VFid_info_file)
    outputfile_1=open(output_file1,"w")
    outputfile_2=open(output_file2,"w")
    with open(blast_filer_file, "r") as fileobject:
        for line in fileobject:
            line=line.strip().split("\t")
            gene_id = line[0]
            contig_name = line[1]
            contig_start = line[2]
            contig_end = line[3]
            strand = line[4]
            sseq = line[11]
            oldline= "\t".join(line[0:11])
            if gene_id == "gene_id":
                #write_header
                print(oldline)
                outputfile_1.write(oldline+"\t"+"species_specific"+"\t"+"plasmid"+"\t"+"ICE"+"\t"+"prophage"+"\t"+"VF_category"+"\n")
            else:
                #read record
                VF_species = info_dict[line[0]][0] 
                VF_plasmid = info_dict[line[0]][1]
                VF_ICE = info_dict[line[0]][2]
                VF_prophage = info_dict[line[0]][3]
                VF_category = info_dict[line[0]][4]
                outputfile_1.write(oldline+"\t"+VF_species+"\t"+VF_plasmid+"\t"+VF_ICE+"\t"+VF_prophage+"\t"+VF_category+"\n")
                index=">"+gene_id+"|"+contig_name +"|"+contig_start+"-"+contig_end + "|" + strand
                outputfile_2.write(index+"\n")
                outputfile_2.write(sseq+"\n")
   
    outputfile_1.close()
    outputfile_2.close()


cal_VFid_draft(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
