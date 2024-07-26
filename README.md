# MetaVF_toolkit

MetaVF_toolkit is a software to profile VF genes in metagenomic sequencing samples

#step1 download MetaVF toolkit

#step2 install MetaVF_toolkit cd /path/to/MetaVF_toolkit

    conda env create -f ./env/MetaVF_tookit.yaml -n MetaVF_toolkit

#step3 decompression database files

    gunzip /path/to/MetaVF_toolkit/databases/pVFGC.fasta.gz
    gunzip /path/to/MetaVF_toolkit/databases/eVFGC.fasta.gz

#step4 add MetaVF_toolkit to environmental variables

    export PATH=/path/to/MetaVF_toolkit:$PATH

#step5 run test

    conda activate MetaVF_toolkit

    metaVF.py -h

    #for draft genemes, please run 

    metaVF.py -p /path/to/MetaVF_toolkit -pjn draft_test -id /path/to/MetaVF_toolkit/example_data/data -o 
    /path/to/MetaVF_toolkit/example_data/test_draft -m draft -c 10 -ti 90 -tc 80

    #attention: the file should ends with ".fna"

    #for pair-end short reads, please run 

    metaVF.py -p /path/to/MetaVF_toolkit -pjn PE_test -id /path/to/MetaVF_toolkit/example_data/data -o 
    /path/to/MetaVF_toolkit/example_data/test_PE -m PE -c 10

    #attention: the file should ends with ".fastq.gz"

    #the output files are located in /path/to/MetaVF_toolkit/example_data/result_test.
