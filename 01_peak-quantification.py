#! /user/bin/python -O
# Jason Matthew Torres
'''
Python 3.5.0 script for merging atac peaks with bedtools
Usage:
module add bedtools
python script
'''

import os,sys,gzip,re
import subprocess as sp
import numpy

work_dir = "/well/mccarthy/users/jason/projects/crds/"
atac_dir = "/well/mccarthy/production/atac-seq/"
peak_dir = atac_dir + "data/human_islets/full_merged_data/peaks/"
bam_dir = atac_dir + "data/human_islets/full_merged_data/bams/"
elife_file = "/well/mccarthy/users/jason/projects/crds/sample_names_100.txt"
out_dir = work_dir + "peak_counts/"
out_file = out_dir + "merged_peaks.txt"
peak_file = "/well/mccarthy/users/maxlouis/CRD_bamfiles/123456_K00000_0123_ABCDEFGHIJ/merged/peaks/R177.rep1-pr.overlap.bfilt.narrowPeak.gz"


def get_sample_list():
    fin = open(elife_file,'r')
    lst = []
    for l in fin:
        lst.append(l.strip())
    fin.close()
    return lst

def pipeline():
    sample_list = get_sample_list()
    print(sample_list)
    fout = open(out_dir+"temp.bed",'w')
    fout.close()
    # Concatenate
    sys.stdout.write("Concatenate files"+"\n")
    # note that ATAC could be replaced with a specific sample name
    command = "zcat " + peak_file + " | awk -v name="+"ATAC"+ " '{print $1,\"\t\",$2,\"\t\",$3,\"\t\",$6,\"\t\",name}'" + " >> " + out_dir + "temp.bed"
    sp.check_call(command,shell=True)
    # Ensure proper formatting
    sys.stdout.write("Ensuring proper formatting"+"\n")
    fout = open(out_dir+"concat_peaks.bed",'w')
    fin = open(out_dir+"temp.bed",'r')
    for line in fin:
        l = line.strip().split()
        fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()
    # Sort
    sys.stdout.write("Sorting"+"\n")
    #command="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools sort -faidx " + atac_dir + "resources/hg19/chr_order_for_merging_peaks.txt -i " + out_dir +"concat_peaks.bed > " + out_dir + "concat_peaks.sorted.bed"
    command="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools sort -faidx " + atac_dir + "resources/hg19/chr_order_for_merging_peaks.txt -i " + out_dir +"concat_peaks.bed > " + out_dir + "concat_peaks.sorted.bed"
    sp.check_call(command,shell=True)
    # Merge (not necessary here)
    #sys.stdout.write("Merging"+"\n")
    #command = "/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools merge -i " + out_dir + "concat_peaks.sorted.bed -c 5 -o collapse > "  + out_dir + "merged_peaks.bed"
    #sp.check_call(command,shell=True)
    # Write SAF file
    sys.stdout.write("Writing SAF file"+"\n")
    fout = open(out_dir+"atac_peaks.saf",'w')
    header = "GeneID\tChr\tStart\tEnd\tStrand\n"
    fout.write(header)
    fin = open(out_dir+"concat_peaks.sorted.bed",'r')
    count=0
    for line in fin:
        count+=1
        l = line.strip().split()
        write_list = ["P"+str(count),l[0],l[1],l[2],"."]
        fout.write("\t".join(write_list)+"\n")
    fin.close()
    fout.close()

    # FeatureCounts
    # "P348"   "HP1507" didn't work at first attempt
    # use instead P328 and HP1507_CMRL
    print("\nFeatureCounts")
    count = 0
    sample_list = ["P328","HP1507_CMRL"]
    for f in sample_list:
        fname = bam_dir + f + ".bam"
        if os.path.isfile(fname) == True:
            count += 1
            sys.stdout.write(f+"\n")
            print count
            command = "/well/mccarthy/production/atac-seq/dependencies/featureCounts -p -a " + out_dir+"atac_peaks.saf" + " -F SAF -o " + out_dir + "atac_peaks.counts-"+f+".txt -T 4 --largestOverlap --ignoreDup -C -f " + fname
            sp.check_call(command,shell=True)


def main():
    pipeline()


if (__name__=="__main__"): main()
