import sys
import os
class VFrecord(object):
    def __init__(self, line):
        fields = line.rstrip().split('\t')
        print(fields)
        self.VF = fields[0]
        self.VFid = fields[1]               
        self.category = fields[3]
        self.species = fields[4]
        self.species_all = fields[5]
        self.length = int(fields[8])
        self.ICE = fields[6]
        self.prophage = fields[7]           
        self.plasmid = fields[9]

class VFcount(object):
    def __init__(self, line):
        fields = line.rstrip().split('\t')
        self.VFid = fields[0]
        self.count = int(fields[1]) 


def cal_VFid(VFid_info_file,rpkm_file,VF_read_count,output_file):
    total_rpk, VF_abun, info_list = cal_rpk(VFid_info_file,VF_read_count,rpkm_file)
       
    outputfile=open(output_file,"w")
    outputfile.write("VFid"+"\t"+"length"+"\t"+"reads_count"+"\t"+"RPK"+"\t"+"TPM"+"\t"+"species_specific"+"\t"+"plasmid"+"\t"+"ICE"+"\t"+"prophage"+"\t"+"VF_category"+"\n")
    for info in info_list:
        VF_TPM=str('%.4g' % float(1000000*VF_abun[info.VFid][1]/total_rpk))
        VF_len=str(info.length)
        VF_reads=str(VF_abun[info.VFid][0])
        VF_rpk=str(round(float(VF_abun[info.VFid][1]),2))
        VF_species =info.species
        VF_category = info.category
        VF_ICE=info.ICE
        VF_prophage=info.prophage
        VF_plasmid=info.plasmid
        VFid= info.VFid
       
        outputfile.write(VFid+"\t"+VF_len+"\t"+VF_reads+"\t"+VF_rpk+"\t"+VF_TPM+"\t"+VF_species+"\t"+VF_plasmid+"\t"+VF_ICE+"\t"+VF_prophage+"\t"+ VF_category +"\n")

    outputfile.close()

def read_VFid_info(VFid_info_file):
    info_list = []
    with open(VFid_info_file,"r") as fileobject:
        next(fileobject)
        for line in fileobject:
            info_list.append(VFrecord(line))
    return info_list

def read_VFid_count(VF_read_count):
    count_list = []
    VFid_detected_list = []
    sum_reads_counts = 0
    with open(VF_read_count,"r") as fileobject:
        for line in fileobject:
            count_list.append(VFcount(line))
            VFid_detected_list.append(VFcount(line).VFid)
            sum_reads_counts += VFcount(line).count
    return count_list,VFid_detected_list,sum_reads_counts

def read_total_counts(rpkm_file):
    with open(rpkm_file,"r") as fileobject:
        for line in fileobject:
            line=line.strip().split("\t")
            if line[0]=="#Reads":
                total_reads=int(line[1])
                
    return total_reads
    
def cal_rpk(VFid_info_file,VF_read_count,rpkm_file):
    info_list = read_VFid_info(VFid_info_file)
    count_list,VFid_detected_list, sum_reads_counts  = read_VFid_count(VF_read_count) 
    total_reads = read_total_counts(rpkm_file)
     
    VF_abun={}
    total_rpk=0
    for VFid_info in info_list:
        if VFid_info.VFid == "Unmapped":
            read_num=total_reads-sum_reads_counts
            rpk=1000*read_num/VFid_info.length
        elif VFid_info.VFid in VFid_detected_list:
            #this has reads count
            read_num = [i.count for i in count_list if i.VFid==VFid_info.VFid][0]
            rpk=1000*read_num/VFid_info.length
        else:
            read_num=0
            rpk=0
        VF_abun.setdefault(VFid_info.VFid,[]).append(read_num)
        VF_abun.setdefault(VFid_info.VFid,[]).append(rpk)
        total_rpk += rpk
    return total_rpk, VF_abun,info_list

    
cal_VFid(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
