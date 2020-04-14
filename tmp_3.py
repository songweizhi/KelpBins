

file = '/Users/songweizhi/Desktop/Kelp_NM/Kelp_NM_dRep99_p32_grouping.txt'


all_identified_taxon = []
kelp_taxon_count_dict = {}
tara_taxon_count_dict = {}
for each in open(file):
    each_split = each.strip().split(',')

    if each_split[2] not in all_identified_taxon:
        all_identified_taxon.append(each_split[2])

    if 'Refined' in each_split[1]:
        if each_split[2] not in kelp_taxon_count_dict:
            kelp_taxon_count_dict[each_split[2]] = 1
        else:
            kelp_taxon_count_dict[each_split[2]] += 1
    else:
        if each_split[2] not in tara_taxon_count_dict:
            tara_taxon_count_dict[each_split[2]] = 1
        else:
            tara_taxon_count_dict[each_split[2]] += 1

print(kelp_taxon_count_dict)
print(tara_taxon_count_dict)
print(all_identified_taxon)


for each_taxon in all_identified_taxon:

    count_in_kelp = 0
    if each_taxon in kelp_taxon_count_dict:
        count_in_kelp = kelp_taxon_count_dict[each_taxon]

    count_in_tara = 0
    if each_taxon in tara_taxon_count_dict:
        count_in_tara = tara_taxon_count_dict[each_taxon]

    print('%s\t%s\t%s' % (each_taxon, count_in_kelp, count_in_tara))








