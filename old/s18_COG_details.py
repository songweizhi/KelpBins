from s0_suspicious_HGTs import suspicious_HGTs


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


# protein_id_cog_file = '/Users/songweizhi/Desktop/GoodBins_0.5_0.05_PG_validated_pcofg_COG_results/protein-id_cog.txt'
# cog_stats_file = '/Users/songweizhi/Desktop/GoodBins_0.5_0.05_PG_validated_pcofg_COG_results/cog_stats.txt'
protein_id_cog_file =   '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_COG_annot/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_protein-id_cog.txt'
cog_stats_file =        '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_COG_annot/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_cog_stats.txt'


suspicious_HGTs_genes = set()
for suspicious_HGT in suspicious_HGTs:
    suspicious_HGT_split = suspicious_HGT.split('___')
    suspicious_HGTs_genes.add(suspicious_HGT_split[0])
    suspicious_HGTs_genes.add(suspicious_HGT_split[1])


cog_cate_to_id_dict = {}
for protein_to_cog in open(protein_id_cog_file):
    protein_id, cog_id, *cog_cates = protein_to_cog.strip().split('\t')
    if protein_id not in suspicious_HGTs_genes:
        for cog_cate in cog_cates:
            if cog_cate not in cog_cate_to_id_dict:
                cog_cate_to_id_dict[cog_cate] = [cog_id]
            else:
                cog_cate_to_id_dict[cog_cate].append(cog_id)


cog_id_to_description_dict = {}
for cog_id in open(cog_stats_file):
    cog_id_split = cog_id.strip().split('\t')
    cog_id_to_description_dict[cog_id_split[0]] = cog_id_split[1]


# cog_cate_to_print = 'E'
# for cog_cate in cog_cate_to_id_dict:
#     if cog_cate == cog_cate_to_print:
#         current_cate_member = cog_cate_to_id_dict[cog_cate]
#         current_cate_member_uniq = unique_list_elements(current_cate_member)
#         for current_cog_id in current_cate_member_uniq:
#             current_cog_id_num = current_cate_member.count(current_cog_id)
#             print('%s\t%s\t%s\t%s' % (cog_cate, current_cog_id, current_cog_id_num, cog_id_to_description_dict[current_cog_id]))


not_print_cog_cate_list = ['E', 'C', 'G', 'P', 'Q']
for cog_cate in cog_cate_to_id_dict:
    if cog_cate not in  not_print_cog_cate_list:
        current_cate_member = cog_cate_to_id_dict[cog_cate]
        current_cate_member_uniq = unique_list_elements(current_cate_member)
        for current_cog_id in current_cate_member_uniq:
            current_cog_id_num = current_cate_member.count(current_cog_id)
            print('%s\t%s\t%s\t%s' % (cog_cate, current_cog_id, current_cog_id_num, cog_id_to_description_dict[current_cog_id]))

