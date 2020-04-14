
BBAY80_81_82_phylum_stats = '/Users/songweizhi/Desktop/Kelp_NM/BBAY80_81_82_phylum.txt'

identified_taxon_list = []
BBAY80_stats_dict = {}
BBAY81_stats_dict = {}
BBAY82_stats_dict = {}
for each in open(BBAY80_81_82_phylum_stats):
    if not each.startswith('Taxon'):

        each_split = each.strip().split('\t')
        phylum_name = each_split[0]
        count_BBAY80 = int(each_split[1])
        count_BBAY81 = int(each_split[2])
        count_BBAY82 = int(each_split[3])

        if phylum_name not in identified_taxon_list:
            identified_taxon_list.append(phylum_name)

        if phylum_name not in BBAY80_stats_dict:
            BBAY80_stats_dict[phylum_name] = count_BBAY80
        else:
            BBAY80_stats_dict[phylum_name] += count_BBAY80

        if phylum_name not in BBAY81_stats_dict:
            BBAY81_stats_dict[phylum_name] = count_BBAY81
        else:
            BBAY81_stats_dict[phylum_name] += count_BBAY81

        if phylum_name not in BBAY82_stats_dict:
            BBAY82_stats_dict[phylum_name] = count_BBAY82
        else:
            BBAY82_stats_dict[phylum_name] += count_BBAY82


print('Taxon\tBBAY80\tBBAY81\tBBAY82\tTotal')
for each_taxon in identified_taxon_list:
    print('%s\t%s\t%s\t%s\t%s' % (each_taxon, BBAY80_stats_dict[each_taxon], BBAY81_stats_dict[each_taxon], BBAY82_stats_dict[each_taxon], BBAY80_stats_dict[each_taxon] + BBAY81_stats_dict[each_taxon] + BBAY82_stats_dict[each_taxon]))




