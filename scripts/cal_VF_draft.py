import sys
#from gene_test import VFrecord
#from gene_test import read_VFid_info

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

def VFid_list(VFid_summary_file):
    VFid_dected=[]
    with open(VFid_summary_file,"r") as fileobject1:
        next(fileobject1)
        for line in fileobject1:
            line=line.strip().split("\t")
            VFid_dected.append(line[0])
    return VFid_dected

def read_VFid_info(VFid_info_file):
    info_list = []
    with open(VFid_info_file,"r") as fileobject:
        next(fileobject)
        for line in fileobject:
            info_list.append(VFrecord(line))
    return info_list

def read_VF_num(VFs_complete_num):
    VF_num={}
    with open(VFs_complete_num,"r") as fileobject2:
        next(fileobject2)
        for line in fileobject2:
            line=line.strip().split("\t")
            VF=line[0]
            complete_num=int(line[1])
            VF_num[VF]=complete_num
    return VF_num


def create_VF_dic(info_list,VFid_dected):
    VF_info_dict={}
    for info in info_list:
        if info.VFid in VFid_dected:
            if info.VF not in VF_info_dict.keys():
                #add value to the dic for the first time
                VF_info_dict.setdefault(info.VF,{})["plasmid"]= list(info.plasmid)
                VF_info_dict.setdefault(info.VF,{})["ICE"]= list(info.ICE)
                VF_info_dict.setdefault(info.VF,{})["prophage"]= list(info.prophage)
                VF_info_dict.setdefault(info.VF,{})["category"]= list(info.category.split("."))
                VF_info_dict.setdefault(info.VF,{})["species_all"]=[s for s in info.species_all.split(";")]
            else:
                VF_info_dict[info.VF]["plasmid"].append(info.plasmid)
                VF_info_dict[info.VF]["ICE"].append(info.ICE)
                VF_info_dict[info.VF]["prophage"].append(info.prophage)
                VF_info_dict[info.VF]["category"].append(info.category)
                VF_info_dict[info.VF]["species_all"] += [s for s in info.species_all.split(";")]

    return VF_info_dict


def filter_redundant(output_file,VF_info_dict,VF_num):
    outputfile=open(output_file,"w")
    outputfile.write("VF"+"\t"+"plasmid"+"\t"+ "ICE"+"\t" + "prophage" +"\t" + "Category" +"\t"+ "Host_species" + "\t" + "Completeness (%)" +"\t" + "Total_genes_in_VF"+"\n" )
    for VF, VF_info in VF_info_dict.items():
        if len(set(VF_info["plasmid"]))==1 and len(set(VF_info["category"]))==1:
            plasmid_new = list(set(VF_info["plasmid"]))[0]
            category_new = list(set(VF_info["category"]))[0]
        else:
            plasmid_new = list(set(VF_info["plasmid"]))[0]
            category_new = list(set(VF_info["category"]))[0]
            print(VF + " background information error!")
            
        if len(set(VF_info["ICE"]))==1:
            ICE_new = list(set(VF_info["ICE"]))[0]
        else:
            ICE_new = "Y"
            
        if len(set(VF_info["prophage"]))==1:
            prophage_new = list(set(VF_info["prophage"]))[0]
        else:
            prophage_new = "Y"
            
        Host_species_new = ';'.join(set(VF_info["species_all"]))
        completeness=round(float(100*len(VF_info["plasmid"]) / int(VF_num[VF])),2)
        Number_of_genes = str(VF_num[VF])
        outputfile.write(VF+"\t"+plasmid_new+"\t"+ ICE_new + "\t" + prophage_new +"\t" + category_new +"\t"+ Host_species_new + "\t" + str(completeness) +"\t" + Number_of_genes +"\n" )

    outputfile.close()


def cal_VFs(VFid_summary_file,VFid_info_file,VFs_complete_num,output_file):

    #step1: fiter VFid_TPM record (filter 0 records)
    VFid_dected = VFid_list(VFid_summary_file)
    #step2 (prepare info file for each VF gene)
    info_list = read_VFid_info(VFid_info_file)
    VF_num = read_VF_num(VFs_complete_num)
    #step3 create VF info dict
    VF_info_dict = create_VF_dic(info_list,VFid_dected)
    #step4 filter and output
    filter_redundant(output_file,VF_info_dict, VF_num)


cal_VFs(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

