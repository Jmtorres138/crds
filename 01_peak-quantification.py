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

work_dir = "/well/mccarthy/users/jason/projects/atac_analyses/evaluate_peaks/"
atac_dir = "/well/mccarthy/production/atac-seq/"
peak_dir = atac_dir + "data/human_islets/full_merged_data/peaks/"
bam_dir = atac_dir + "data/human_islets/full_merged_data/bams/"
elife_file = "/well/mccarthy/users/jason/projects/atac_analyses/elife_samples_names_for_jason.txt"
out_dir = work_dir + "eLife2018/"
out_file = out_dir + "merged_peaks.txt"
suffix = ".tn5.pval0.01.300K.bfilt.narrowPeak.gz"

def get_elife_list():
    fin = open(elife_file,'r')
    lst = []
    for l in fin:
        lst.append(l.strip())
    fin.close()
    return lst

def pipeline():
    elife_list = get_elife_list()
    elife_list.append("HP1535")
    elife_list.append("HP1507_CMRL")
    print(elife_list)
    fout = open(out_dir+"temp.bed",'w')
    fout.close()
    # Concatenate
    for f in elife_list:
        fname = peak_dir + f + suffix
        if os.path.isfile(fname) == True:
            sys.stdout.write(f+"\n")
            #command = "zcat " + fname + " | cut -f 1,2,3,6 >> " + out_dir + "concat_peaks.bed"
            command = "zcat " + fname + " | awk -v name="+f+ " '{print $1,\"\t\",$2,\"\t\",$3,\"\t\",$6,\"\t\",name}'" + " >> " + out_dir + "temp.bed"
            sp.check_call(command,shell=True)
    # Ensure proper formatting
    fout = open(out_dir+"concat_peaks.bed",'w')
    fin = open(out_dir+"temp.bed",'r')
    for line in fin:
        l = line.strip().split()
        fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()
    # Sort
    command="/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools sort -faidx " + atac_dir + "resources/hg19/chr_order_for_merging_peaks.txt -i " + out_dir +"concat_peaks.bed > " + out_dir + "concat_peaks.sorted.bed"
    #command = "sort -k1,1 -k2,2n " +  out_dir + "concat_peaks.bed" + " | uniq -u " " > " + out_dir + "concat_peaks.sorted.bed"
    sp.check_call(command,shell=True)
    # Merge
    command = "/apps/well/bedtools/2.24.0-18-gb0bc5b7/bin/bedtools merge -i " + out_dir + "concat_peaks.sorted.bed -c 5 -o collapse > "  + out_dir + "merged_peaks.bed"
    sp.check_call(command,shell=True)
    # Write SAF file
    fout = open(out_dir+"merged_peaks.saf",'w')
    header = "GeneID\tChr\tStart\tEnd\tStrand\n"
    fout.write(header)
    fin = open(out_dir+"merged_peaks.bed",'r')
    count=0
    for line in fin:
        count+=1
        l = line.strip().split()
        write_list = ["P"+str(count),l[0],l[1],l[2],"."]
        fout.write("\t".join(write_list)+"\n")
    fin.close()
    fout.close()

    # FeatureCounts
    print("\nFeatureCounts")
    for f in elife_list:
        fname = bam_dir + f + ".bam"
        if os.path.isfile(fname) == True:
            sys.stdout.write(f+"\n")
            command = "/well/mccarthy/production/atac-seq/dependencies/featureCounts -p -a " + out_dir+"merged_peaks.saf" + " -F SAF -o " + out_dir + "merged_peaks.counts-"+f+".txt -T 4 --largestOverlap --ignoreDup -C -f " + fname
            sp.check_call(command,shell=True)



def main():
    pipeline()
    #merge_bed_file(out_file)


if (__name__=="__main__"): main()
