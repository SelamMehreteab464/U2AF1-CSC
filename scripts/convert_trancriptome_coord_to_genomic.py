import pandas as pd
import numpy as np


transcriptomebed = "/private/groups/brookslab/smehrete/RNA_modification/FLAIR/07082026_flair_combined_transcriptome/07082026_flair_combined_transcriptome.bed"
myfile = "/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/m6anetoutputs/data.site_proba.csv"
output = "/private/groups/brookslab/smehrete/RNA_modification/07102026.m6anet.results/m6anetoutputs/m6anet_data.site_proba.genomic.csv"


isotoblocks = {}
geneloc= {}
for line in open(transcriptomebed):
    line = line.rstrip().split('\t')
    thischr, iso, dir, start, esizes, estarts = line[0], line[3], line[5], int(line[1]), [int(x) for x in line[10].split(',')[:-1]], [int(x) for x in line[11].split(',')[:-1]]
    exonblocks = []  ##block is gstart, tstart, len
    if dir == '+':
        currtstart = 0
        for i in range(len(esizes)):
            exonblocks.append((currtstart, estarts[i] + start, esizes[i]))
            currtstart += esizes[i]
    else:
        currtstart = sum(esizes)
        for i in range(len(esizes)):
            exonblocks.append((currtstart, estarts[i] + start, esizes[i]))
            currtstart -= esizes[i]
    isotoblocks[iso] = (thischr, dir, exonblocks)

def convert_pos_to_genomic(isoinfo, tpos):
    thischr, dir, blocks = isoinfo #isotoblocks[iso]
    if dir == '-':
        tpos -= 1
    genomePos = None
    if dir == '+':
        for tstart, gstart, bsize in blocks:
            if tstart <= tpos < tstart + bsize:
                genomePos = gstart + (tpos - tstart)
    else:
        for tstart, gstart, bsize in blocks:
            if tstart - bsize <= tpos < tstart:
                genomePos = gstart + (tstart - tpos)
    genomePos -= 1
    return genomePos, thischr

with open(output, 'w') as out:
    out.write('transcript_id\ttranscript_position\tread_index\tprobability_modified\tchrom\tgenomepos\n')

    with open(myfile) as f:
        next(f)
        for line in f:
            line = line.rstrip().split(',')
            iso, tpos = line[0], int(line[1]) + 1
            genomepos, chrom = convert_pos_to_genomic(isotoblocks[iso], tpos)
            #with open(output, 'a') as out:
            out.write('\t'.join(line) + f'\t{chrom}\t{genomepos}\n')

print('converted')

   ##DO STUFF, OUTPUT TO NEW FILE



