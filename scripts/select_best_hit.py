import sys
class BlastHit(object):
    def __init__(self, line):
        fields = line.rstrip().split('\t')
        self.gene_id = fields[0]                #query_id
        self.contig_name = fields[1]            #subject_id
        self.contig_start = int(fields[20])     #subject_start
        self.contig_end = int(fields[21])       #subject_end
        self.gene_id_long = fields[4]           #query_id_full
        self.pcid = float(fields[5])            #identity in percentage
        self.alignment_length = int(fields[6])  #length of alignment
        self.score = float(fields[12])          #bit score
        self.qcovs = float(fields[14])          #query coverage
        self.typep = fields[18]                 #type of alleles: perfect or strict
        self.sseq_new = fields[19]              #subject sequence
        self.strand = fields[16]                #subject strand


def select_best_blast_hit(blast_file, output_file):
    blast_hits = [] #all_records
    with open(blast_file,"r") as fileobject:
        next(fileobject)
        for line in fileobject:
            blast_hits.append(BlastHit(line))


    # Clean up redundant hits so we only have one hit for each part of the genome.
    blast_hits = cull_redundant_hits(blast_hits)
    #return blast_hits
    f = open(output_file,"w")
    f.write(str("gene_id")+"\t"+\
                "contig_name" + "\t" + \
                "contig_start"+"\t"+\
                "contig_end" + "\t" +\
                "strand" + "\t" +\
                "gene_id_long" + "\t"+\
                "pcid" + "\t"+\
                "qcovs"+ "\t"+\
                "score" + "\t"+\
                "alignment_length" + "\t"+\
                "typep" + "\t"+\
                "sseq_new" + "\t"+\
                "\n")
    for blast_hit in blast_hits:
        f.write(str(blast_hit.gene_id)+"\t"+\
                blast_hit.contig_name + "\t" + \
                str(blast_hit.contig_start)+"\t"+\
                str(blast_hit.contig_end) + "\t" +\
                str(blast_hit.strand) + "\t" +\
                str(blast_hit.gene_id_long) + "\t"+\
                str(blast_hit.pcid) + "\t"+\
                str(blast_hit.qcovs) + "\t"+\
                str(blast_hit.score) + "\t"+\
                str(blast_hit.alignment_length) + "\t"+\
                str(blast_hit.typep) + "\t"+\
                str(blast_hit.sseq_new) + "\t"+\
                "\n")


    f.close()



def cull_redundant_hits(blast_hits):
    """
    Cull out redundant hits
    """
    # Sort the hits from best to worst. Hit quality is defined as the product of gene coverage,
    # identity and score.
    blast_hits = sorted(blast_hits, key=lambda x: (1/(x.pcid * x.score * x.qcovs), x.gene_id))

    filtered_blast_hits = []

    for h in blast_hits:
        if not overlapping(h, filtered_blast_hits):
            filtered_blast_hits.append(h)

    return filtered_blast_hits


def overlapping(hit, existing_hits):
    # Only consider hits in the same reading frame.
    existing_hits = [h for h in existing_hits if
                     h.strand == hit.strand and
                     h.contig_name == hit.contig_name]

    for existing_hit in existing_hits:
        if hits_overlap(hit, existing_hit):
            return True

    return False


def hits_overlap(a, b):
    if a.contig_start <= b.contig_end and b.contig_start <= a.contig_end:  # There is some overlap
        allowed_overlap = 10
        overlap_size = len(range(max(a.contig_start, b.contig_start),
                                 min(a.contig_end, b.contig_end) + 1))
        return overlap_size > allowed_overlap
    else:
        return False


select_best_blast_hit(sys.argv[1],sys.argv[2])

