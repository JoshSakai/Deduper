#!/usr/bin/env python

import re
import argparse

def get_args():
    """argparse input files"""
    parser = argparse.ArgumentParser(description='file name to be used')
    parser.add_argument('-f', '--file', type=str, required=True, help='File path to sorted sam file')
    parser.add_argument('-p', '--paired', required=False, default=False, help='designates file is paired end (not single-end)')
    parser.add_argument('-u', '--umi', type=str, required=False, default=None, help= 'designates file containing the list of UMIs (unset if randomers instead of UMIs)')
    #parser.add_argument('-h', '--help', required=False, help='prints a help message')
    return parser.parse_args()

args =get_args()
file = args.file
paired = args.paired
umi_flag = args.umi

def umi_ref(umi):
    """given the input file from the argparse this takes the umis within it and creates a reference array as a key for deduplicating"""
    """importantly the input file must have am umi on each line"""

    with open(umi,'r') as u_ref :
        umi_list = []
        for line in u_ref:
            y = line.strip().split()[0]
            umi_list.append(y)
        return umi_list

def deduper(file):
    """
key factors that need to account for : the adjusted 5' postion, the sequence, and the cigar string, the bitwise flag, and the QNAME for soft clipping

First
    open up the list of umis and store them in an array
    open both sam files and store a group of reads by chromosome or scaffold to prevent too many reads from being stored at any one time.
    open up 3 outfiles one file for UMI mismatches, one for PCR duplicates, and one for the rest of the good reads

Second
    check the reads UMI's to see if they match the reference, if they don't write them to a bad out file
    iterate through the stored reads checking the 5' most position accounting for soft clipping by subtracting the given value by S if the cigar string indicates there was left-most softclipping.

Third
    what is the alignment of the new read (+/-)?
    Any reads that have D in the cigar string AND are (-) must have their starting position adjusted similar to accounting for right-soft clipping
    Deletions and insertion do not make an impact on (+) sequence positions

Fourth
    if a read at the leftmost position has not been recorded, write it to an outfile then move on to the next stored read.

Fifth
    if the next/subsequent read has the same left-most position then check to see if the UMI's attached to the reads are the same

        if yes then toss the read with the lower mean quality score since these are PCR duplicates
        if no then it's not a PCR duplicate, and it can be written to the outfile

    if the next sequence starts at a different position than the last observed write it to the outfile and move on to the next read
"""
    #making a set to store the references
    dupe_check = set()

    #making a counter of the number of duplicates
    dupe_count = 0

    #setting the last rname to none for the start of the read
    last_rname = None

    with open(file,"r") as f:          #open the input sam file
        with open(file+'_deduped.sam',"w") as out:         #create an output sam file for the good reads
            for line in f:
            #sam files have header lines so we need to output those into our deduped file first
                if line.startswith('@'):
                    out.write(line)
                    continue
                    #moving on to looking at actual aligned reads
                read = line.strip()
                    #grabbing all relavent information from the line reads
                qname = read.split('\t')[0]
                this_umi = qname.split(':')[-1] # the umi is the last thing in the qname
                flag = int(read.split('\t')[1])
                rname = read.split('\t')[2]
                pos = int(read.split('\t')[3])
                #qual = read.split('\t')[10]
                cig = read.split('\t')[5]

                if this_umi in u_ref:
                #taking in the initial position given from the headerline and the cigar string this function adjusts the position of reads being compared
                #for either right or left soft clipping this returns the 'real position' of the 5' most position of the read
                    #start with initial blank values for left or right clipping
                    lclip = 0
                    rclip = 0

                    cig_fw = re.findall(r'(\d+)([M,N,S]{1})', cig)
                    cig_rv = re.findall(r'(\d+)([M,N,S,D]{1})', cig)

                    fw_lengths = [int(item[0]) for item in cig_fw]
                    fw_len = sum(fw_lengths)

                    rv_lengths = [int(item[0]) for item in cig_rv]
                    rv_len = sum(rv_lengths)

                    for cha in cig_rv:
                        if cha[-1] == 'S':
                            rclip = int(cig_rv[-1][0])
                        if cig_rv[0][1] == 'S':
                            lclip =  int(cig_rv[0][0])

                    #for each read in the group files
                    rev_comp = False
                    if ((flag & 16) == 16):
                        rev_comp = True
                        #set a flag to see which strand the read is on
                    if rev_comp is False:  #if this read is on the (+) strand
                        real_start = pos - lclip   # adjusting for left clipping
                        end = real_start + fw_len   # finding the end position given length
                        real_end = end + rclip   # adjusting for right clipping

                        if (this_umi, real_start, rev_comp) in dupe_check:
                            dupe_count += 1
                        else:
                            dupe_check.add((this_umi, real_start, rev_comp))
                            out.write(line)

                    else:   # if this read is on the (-) strand
                        real_end = pos - lclip   # adjusting for left clipping
                        start = real_end + rv_len # finding the end position given length
                        real_start = start + rclip   # adjusting for right clipping

                        if (this_umi, real_start, rev_comp) in dupe_check:
                                dupe_count += 1
                        else:
                            dupe_check.add((this_umi, real_start, rev_comp))
                            out.write(line)

                if last_rname is None:
                    last_rname = rname
                elif last_rname == rname:
                    pass
                else:
                    print(str(dupe_count)+' PCR duplicates removed from '+str(last_rname))
                    dupe_check.clear()
                    dupe_count = 0
                    last_rname = rname
                continue
                
    print(str(dupe_count)+' PCR duplicates removed from '+str(last_rname))

"""Now we have defined all the functions we can actually run our code"""

#for testing purposes

if umi_flag is not None:         #if we get an actual umi file run the umi ref function to make a list of umis
    u_ref = umi_ref(umi_flag)
else:
    u_ref=[]
    #info = file.readlines().split('\t')[0]
    #u_ref = file.readlines().split(':')[-1]
    with open(file,"r") as f:
            for line in f:
            #sam files have header lines so we need to output those into our deduped file first
                if line.startswith('@'):
                    continue
                    #moving on to looking at actual aligned reads
                read = line.strip()
                    #grabbing all relavent information from the line reads
                qname = read.split('\t')[0]
                this_umi = qname.split(':')[-1] # the umi is the last thing in the qname
                u_ref.append(this_umi)

deduper(file)