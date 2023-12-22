#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import argparse
import subprocess
import re
from functools import reduce

def write_yaml(dic_file_need_profile_direct_path, yaml_file_name, project_name, outdir_direct_path, mode, ti, tc):
    all_file_names = os.listdir(dic_file_need_profile_direct_path)
    list_names = [re.match('(.*).fna', name).groups()[0] for name in all_file_names if mode=="draft" and re.match('(.*).fna$', name)]
    if len(list_names) == 0:
        list_names = [re.match('(.*)_[1,2].fastq.gz', name).groups()[0] for name in all_file_names if mode=="PE" and re.match('(.*)_[1,2].fastq.gz$', name)]
    yaml_file = open(yaml_file_name, "w")
    yaml_lines = "#snakemake\nmetaVF:\n    sample: %s\n    project: %s\n    datadir: %s/\n    outdir: %s\n    thre_identity: %d\n    thre_coverage: %d"%\
            ("[%s]"%(reduce(lambda x,y: "%s,%s"%(x,y), list_names))\
            , project_name, dic_file_need_profile_direct_path, outdir_direct_path, ti, tc)
    yaml_file.write(yaml_lines)

def main():
    parser = argparse.ArgumentParser(description='Profile VF genes in the metagenomic sequencing sample')

    parser.add_argument('-p',type=str, dest='path', action='store', required = False,
                        default=".",
                        help='direct path of MetaVF_toolkit (example: /data/home/MetaVF_toolkit)')
    parser.add_argument('-pjn','--project_name', type=str, dest='project_name', action='store', required = False,
                        default="metavfjob",
                        help='name of the project(example: draft)')
    parser.add_argument('-id','--in_data_dir', type=str, dest='direct_indata_path', required=True, \
            help='the direct path of the in data directory(example: /data/home/MetaVF_toolkit/example/data/draft)')
    parser.add_argument('-outd','--output_dir', type=str, dest='direct_output_path', required=False, \
            default='', help='the direct path of the output files directory(example: /data/home/MetaVF_toolkit/example/data)')
    parser.add_argument('-ti','--thre_identity', type=int, required=False, default=90, \
            help='identity to filter blastn hits, in draft mode only') 
    parser.add_argument('-tc', '--thre_coverage', type=int, required=False, default=80, \
            help='coverage to filter blastn hits, in draft mode only') 
    parser.add_argument('-m','--mode', type = str, dest='mode', action='store', required = True,
                        choices = ["PE","draft"],
                        help='type of input data (choose from "PE","draft")')
    parser.add_argument('-c', '--cores',type = int, dest='cores', action='store', required = False,
                        default="20",
                        help='the max core of snakemake (default: 20)')
    args = parser.parse_args()

    # for the example id, father_dir_of_the_project is /data/home/MetaVF_toolkit/example/data
    father_dir_of_the_project = args.direct_indata_path
    if args.direct_output_path == "":
        direct_output_path = "%s/%s_result/"%(father_dir_of_the_project, args.project_name)
    elif args.direct_output_path != "":
        direct_output_path = args.direct_output_path
    yaml_file_name = "%s/%s_config.yaml"%(direct_output_path, args.project_name)
    
    if not os.path.exists('%s/'%direct_output_path):
        os.mkdir('%s/'%direct_output_path)
    write_yaml(args.direct_indata_path, yaml_file_name, args.project_name, direct_output_path, args.mode, args.thre_identity, args.thre_coverage)

    os.chdir(args.path)
    if args.mode=="SE":
        os.system("snakemake -s ./workflows/Snakefile_SE --configfile %s --use-conda --rerun-incomplete --cores  %d -p" %(yaml_file_name, args.cores))
    elif args.mode=="PE":
        os.system("snakemake -s ./workflows/Snakefile_PE --configfile %s --use-conda --rerun-incomplete --cores %d -p" %(yaml_file_name, args.cores))
    elif args.mode=="draft":
        os.system("snakemake -s ./workflows/Snakefile_draft --configfile %s --use-singularity --rerun-incomplete --cores %d -p" %(yaml_file_name, args.cores))
    os.remove(yaml_file_name)

if __name__ == "__main__":
    main()
