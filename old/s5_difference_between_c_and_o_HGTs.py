
HGT_PG_validated_c = '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_c15_HGTs_PG_validated.txt'
HGT_PG_validated_o = '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_o34_HGTs_PG_validated.txt'


hgt_o_gene_list = set()
hgt_o_gene1_list = set()
donor_o_list = set()
recipient_o_list = set()
for hgt_o in open(HGT_PG_validated_o):
    if not hgt_o.startswith('Gene_1'):
        hgt_o_split = hgt_o.strip().split('\t')
        donor_o_genome = hgt_o_split[7].split('-->')[0]
        recipient_o_genome = hgt_o_split[7].split('-->')[1]

        donor_o_gene = hgt_o_split[0]
        recipient_o_gene = hgt_o_split[1]
        if donor_o_genome in hgt_o_split[1]:
            donor_o_gene = hgt_o_split[1]
            recipient_o_gene = hgt_o_split[0]

        hgt_o_gene_list.add(hgt_o_split[0])
        hgt_o_gene1_list.add(hgt_o_split[0])
        hgt_o_gene_list.add(hgt_o_split[1])
        donor_o_list.add(donor_o_gene)
        recipient_o_list.add(recipient_o_gene)


hgt_c_gene_list = set()
hgt_c_gene1_list = set()
donor_c_list = set()
recipient_c_list = set()
for hgt_c in open(HGT_PG_validated_c):
    if not hgt_c.startswith('Gene_1'):
        hgt_c_split = hgt_c.strip().split('\t')
        donor_c_genome = hgt_c_split[7].split('-->')[0]
        recipient_c_genome = hgt_c_split[7].split('-->')[1]

        donor_c_gene = hgt_c_split[0]
        recipient_c_gene = hgt_c_split[1]
        if donor_c_genome in hgt_c_split[1]:
            donor_c_gene = hgt_c_split[1]
            recipient_c_gene = hgt_c_split[0]

        hgt_c_gene_list.add(hgt_c_split[0])
        hgt_c_gene1_list.add(hgt_c_split[0])
        hgt_c_gene_list.add(hgt_c_split[1])
        donor_c_list.add(donor_c_gene)
        recipient_c_list.add(recipient_c_gene)


shared = 0
for each in hgt_c_gene_list:
    if each in hgt_o_gene_list:
        shared += 1


shared_gene1 = 0
for each in hgt_c_gene1_list:
    if each in hgt_o_gene1_list:
        shared_gene1 += 1


shared_d = 0
for each in donor_c_list:
    if each in donor_o_list:
        shared_d += 1

shared_r = 0
for each in recipient_c_list:
    if each in recipient_o_list:
        shared_r += 1

print('\to\tboth\tc')
print('all\t%s\t%s\t%s' % (len(hgt_o_gene_list), shared, len(hgt_c_gene_list)))
print('g1\t%s\t%s\t%s' % (len(hgt_o_gene1_list), shared_gene1, len(hgt_c_gene1_list)))

#print('d\t%s\t%s\t%s' % (len(donor_o_list), shared_d, len(donor_c_list)))
#print('r\t%s\t%s\t%s' % (len(recipient_o_list), shared_r, len(recipient_c_list)))





