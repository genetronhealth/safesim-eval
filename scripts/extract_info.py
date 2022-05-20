#!/usr/bin/env python

import sys, pysam

vcf = sys.argv[1]
dataclass = sys.argv[2]

bed_list = []

if dataclass == 'WES':
    gs = sys.path[0] + '/../simfiles/wes_simmut_gs.txt'
    bed = sys.path[0] + '/../simfiles/wes_simmut.bed'
else:
    gs = sys.path[0] + '/../simfiles/ctDNA_simmut_gs.txt'
    bed = sys.path[0] + '/../simfiles/' + dataclass + '_simmut.bed'

def print_info(c, mut, af, dp):
    if c == 'WES':
        print ('{}\t{}\t{}'.format(mut, af, dp))
    else:
        print ('{}\t{}'.format(mut, af))
        
with open(bed, 'r') as b:
    for eachline in b:
        bed_list.append(eachline.strip())

vcffile = pysam.VariantFile(vcf)
with open(gs, 'r') as f:
    for eachline in f:
        af, dp = 0, 0
        inbed = False
        mut, gs_af = eachline.strip().split('\t')[0:2]
        chr, pos, alt = mut.split('-')
        for seg in bed_list:
            seg_c, seg_s, seg_e = seg.split('\t')
            if chr == seg_c and int(seg_s) < int(pos) <= int(seg_e):
                inbed = True
        if not inbed:
            af, dp = 'NA', 'NA'
            print_info(dataclass, mut, af, dp)
            continue
        if not list(vcffile.fetch(region = chr + ':' + pos + '-' + pos)):
            print_info(dataclass, mut, af, dp)
            continue
        is_info_printed = False
        maxDP = 0
        for rec in vcffile.fetch(region = chr + ':' + pos + '-' + pos):
            if 'DP' in rec.samples[0]: maxDP = max([maxDP, rec.samples[0]['DP']])
            if rec.alts[0] == alt:
                ad = rec.samples[0]['AD'][1]
                dp = rec.samples[0]['DP']
                af = float(ad)/dp
                print_info(dataclass, mut, af, dp)
                is_info_printed = True
                break
        if not is_info_printed: 
            print_info(dataclass, mut, af, dp)
            #print_info(dataclass, mut, 0, maxDP)
     
