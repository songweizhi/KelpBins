
def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


cog_annotation_results_wd = '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_COG_annot'
pwd_cog_stats_file =        '%s/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_cog_stats.txt'      % cog_annotation_results_wd
pwd_func_stats_file =       '%s/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_func_stats.txt'     % cog_annotation_results_wd
pwd_protein_id_cog_file =   '%s/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_protein-id_cog.txt' % cog_annotation_results_wd
interested_COG_category = 'G'


# get COG id to description dict
COG_id_to_description_dict = {}
for each_id in open(pwd_cog_stats_file):
    each_id_split = each_id.strip().split('\t')
    COG_id_to_description_dict[each_id_split[0]] = each_id_split[1]


interested_COG_id_list = []
for each_gene in open(pwd_protein_id_cog_file):
    each_gene_split = each_gene.strip().split()
    each_gene_COGs = each_gene_split[2:]
    if interested_COG_category in each_gene_COGs:
        interested_COG_id_list.append(each_gene_split[1])


interested_COG_id_list_uniq = unique_list_elements(interested_COG_id_list)

for each_uniq in interested_COG_id_list_uniq:
    for_print = '%s\t%s\t%s\t%s' % (each_uniq, interested_COG_id_list.count(each_uniq), float("{0:.2f}".format(interested_COG_id_list.count(each_uniq)*100/len(interested_COG_id_list))), COG_id_to_description_dict[each_uniq])
    print(for_print)

