import pysam


folder = '/private/groups/brookslab/smehrete/filtered_align_bam/'
files = [
    'M1C2_S29.hg38aligned.filtered.bam',
    'M1C3_S37.hg38aligned.filtered.bam',
    'M1D2_S28.hg38aligned.filtered.bam',
    'M1D3_S36.hg38aligned.filtered.bam',
    'M2C2_S27.hg38aligned.filtered.bam',
    'M2C3_S35.hg38aligned.filtered.bam',
    'M2D2_S26.hg38aligned.filtered.bam',
    'M2D3_S34.hg38aligned.filtered.bam',
    'W1C2_S25.hg38aligned.filtered.bam',
    'W1C3_S33.hg38aligned.filtered.bam',
    'W1D2_S24.hg38aligned.filtered.bam',
    'W1D3_S32.hg38aligned.filtered.bam',
    'W2C2_S23.hg38aligned.filtered.bam',
    'W2C3_S31.hg38aligned.filtered.bam',
    'W2D2_S22.hg38aligned.filtered.bam',
    'W2D3_S30.hg38aligned.filtered.bam'
]

def get_introns(a):
    ref = a.reference_start
    introns = []
    for block in a.cigartuples:
        if block[0] in {0, 7, 8}:  # match, consumes both
            ref += block[1]
        elif block[0] in {2, 3}:  # consumes reference #2 is deletion
            if block[0] == 3:  # intron
                introns.append((ref, ref + block[1]))
            ref += block[1]
    return introns

#dsa, usa = 31722330, 31722348
#e1A, e1Bds, e1Bus, e1C = 31722617, 31723734, 31723863, 31722416
intron_to_counts = {}
for file in files:
#    e1A_c, e1Bds_c, e1Bus_c, e1C_c = 0, 0, 0, 0
    for a in pysam.AlignmentFile(folder + file, 'rb').fetch('20', 31721346, 31724156):
        if a.is_unmapped or a.cigartuples is None:
            continue
        introns = get_introns(a)
        if len(introns) >= 1:
            first_intron = introns[-1]
            if first_intron not in intron_to_counts:
                intron_to_counts[first_intron] = 0
            intron_to_counts[first_intron] += 1
for intron, count in intron_to_counts.items():
    if count >= 10:
        print(f'chr20:{intron[0]}-{intron[1]}', count)
    #my_acceptor, my_donor = first_intron
            #if my_acceptor == dsa or my_acceptor == usa:
            # if my_acceptor == dsa:
            #if my_acceptor == usa:
             #   if my_donor == e1A:
              #      e1A_c += 1
               # if my_donor == e1Bds:
                #    e1Bds_c += 1
                #if my_donor == e1Bus:
                 #   e1Bus_c += 1
                #if my_donor == e1C:
                 #   e1C_c += 1
    #print(file.split('.')[0], e1A_c, e1Bds_c, e1Bus_c, e1C_c)


# intron_to_name = {
#     (31722330, 31722617):'e1A-dsa',
#     (31722348, 31722617):'e1A-usa',
#     (31722330, 31723734):'e1Bds-dsa',
#     (31722330, 31723863):'e1Bus-dsa',
#     (31722348, 31723734):'e1Bds-usa',
#     (31722348, 31723863):'e1Bus-usa',
#     (31722330, 31722416):'e1C-dsa',
# }

# """
# (31722330, 31722617) 955
# (31722348, 31722617) 1017
# (31722330, 31723734) 295
# (31722330, 31723847) 16
# (31722330, 31723863) 279
# (31722348, 31723734) 142
# (31722330, 31722443) 36
# (31722348, 31723863) 110
# (31689423, 31721654) 5
# (31722330, 31722416) 117
# (31666086, 31721654) 10
# (31722330, 31722476) 9
# (31722348, 31722443) 8
# """
# for intron, count in intron_to_counts.items():
#     if count >= 5:
#         print(f'chr20:{intron[0]}-{intron[1]}', count)
