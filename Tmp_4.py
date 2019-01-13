import os

wd = '/Users/songweizhi/Desktop/KelpBins/difference_between_c_and_o_HGTs'

os.chdir(wd)

#HGT_PG_validated_c =        'GoodBins_0.5_0.05_c15_HGTs_PG_validated.txt'
#HGT_PG_validated_o =        'GoodBins_0.5_0.05_o34_HGTs_PG_validated.txt'

HGT_PG_validated_c =        'GoodBins_0.5_0.05_c15_HGTs_BM.txt'
HGT_PG_validated_o =        'GoodBins_0.5_0.05_o34_HGTs_BM.txt'
taxon_to_group_id_file_c =  'GoodBins_0.5_0.05_c15_group_to_taxon.txt'
taxon_to_group_id_file_o =  'GoodBins_0.5_0.05_o34_group_to_taxon.txt'
gtdbtk_output =             'gtdbtk.bac120.classification_528_r86.tsv'


blast_hit_in_one_line_c =                       'GoodBins_0.5_0.05_c15_3_blastn_results_filtered_in_one_line.tab'
blast_hit_in_one_line_o =                       'GoodBins_0.5_0.05_o34_3_blastn_results_filtered_in_one_line.tab'
blast_hit_in_one_line_extracted =               'GoodBins_0.5_0.05_o34_3_blastn_results_filtered_in_one_line_extracted.tab'
blast_hit_in_one_line_extracted_tmp1 =          'GoodBins_0.5_0.05_o34_3_blastn_results_filtered_in_one_line_extracted_tmp1.tab'
blast_hit_in_one_line_extracted_tmp1_sorted =   'GoodBins_0.5_0.05_o34_3_blastn_results_filtered_in_one_line_extracted_tmp1_sorted.tab'
blast_hit_in_one_line_extracted_tmp1_sorted2 =  'GoodBins_0.5_0.05_o34_3_blastn_results_filtered_in_one_line_extracted_tmp1_sorted_renamed.tab'


hgt_o_gene1_list = set()
for hgt_o in open(HGT_PG_validated_o):
    if not hgt_o.startswith('Gene_1'):
        hgt_o_split = hgt_o.strip().split('\t')
        hgt_o_gene1_list.add(hgt_o_split[0])

hgt_c_gene1_list = set()
for hgt_c in open(HGT_PG_validated_c):
    if not hgt_c.startswith('Gene_1'):
        hgt_c_split = hgt_c.strip().split('\t')
        hgt_c_gene1_list.add(hgt_c_split[0])


shared_gene1 = set(hgt_c_gene1_list).intersection(hgt_o_gene1_list)


hgt_c_gene1_list_uniq = []
for each in hgt_c_gene1_list:
    if each not in shared_gene1:
        hgt_c_gene1_list_uniq.append(each)

hgt_o_gene1_list_uniq = []
for each in hgt_o_gene1_list:
    if each not in shared_gene1:
        hgt_o_gene1_list_uniq.append(each)

print('uniq_c\tshare\tuniq_o')
print('%s\t%s\t%s' % (len(hgt_c_gene1_list_uniq), len(shared_gene1), len(hgt_o_gene1_list_uniq)))


# blast_hit_in_one_line_extracted_handle = open(blast_hit_in_one_line_extracted, 'w')
# for each_query_c in open(blast_hit_in_one_line_c):
#     each_query_c_split = each_query_c.strip().split('\t')
#     each_query_c_gene_id = each_query_c_split[0].split('|')[1]
#     if each_query_c_gene_id in hgt_c_gene1_list_uniq:
#         blast_hit_in_one_line_extracted_handle.write('uniq_c_in_c\t%s' % each_query_c)
#     if each_query_c_gene_id in shared_gene1:
#         blast_hit_in_one_line_extracted_handle.write('both_2_in_c\t%s' % each_query_c)
#     if each_query_c_gene_id in hgt_o_gene1_list_uniq:
#         blast_hit_in_one_line_extracted_handle.write('uniq_o_in_c\t%s' % each_query_c)
#
# for each_query_o in open(blast_hit_in_one_line_o):
#     each_query_o_split = each_query_o.strip().split('\t')
#     each_query_o_gene_id = each_query_o_split[0].split('|')[1]
#     if each_query_o_gene_id in hgt_o_gene1_list_uniq:
#         blast_hit_in_one_line_extracted_handle.write('uniq_o_in_o\t%s' % each_query_o)
#     if each_query_o_gene_id in shared_gene1:
#         blast_hit_in_one_line_extracted_handle.write('both_2_in_o\t%s' % each_query_o)
#     if each_query_o_gene_id in hgt_c_gene1_list_uniq:
#         blast_hit_in_one_line_extracted_handle.write('uniq_c_in_o\t%s' % each_query_o)
#
# blast_hit_in_one_line_extracted_handle.close()



blast_hit_in_one_line_extracted_tmp1_handle = open(blast_hit_in_one_line_extracted_tmp1, 'w')
for each_extracted in open(blast_hit_in_one_line_extracted):
    each_extracted_split = each_extracted.strip().split()

    query_split = each_extracted_split[1].split('|')
    needed = [query_split[1], each_extracted_split[0], query_split[0].split('_')[0]]
    for each_hit in each_extracted_split[2:]:
        each_hit_split = each_hit.split('|')
        each_hit_needed = '%s|%s' % (each_hit_split[0].split('_')[0], each_hit_split[2])
        needed.append(each_hit_needed)

    # get group list for each query
    query_hits_group_iden_dict = {}
    for each_hit in needed[3:]:
        each_hit_group = each_hit.split('|')[0]
        each_hit_iden = float(each_hit.split('|')[1])
        if each_hit_group not in query_hits_group_iden_dict:
            query_hits_group_iden_dict[each_hit_group] = [each_hit_iden]
        else:
            query_hits_group_iden_dict[each_hit_group].append(each_hit_iden)

    query_hits_group_iden_dict_mean = {}
    for each_key in query_hits_group_iden_dict:
        iden_list = query_hits_group_iden_dict[each_key]
        iden_mean = sum(iden_list)/len(iden_list)
        #query_hits_group_iden_dict_mean[each_key] = '%s' % float("{0:.3f}".format(iden_mean))
        query_hits_group_iden_dict_mean[each_key] = '%s(%s)' % (float("{0:.3f}".format(iden_mean)), len(iden_list))

    self_group = query_split[0].split('_')[0]
    self_group_mean_iden = 'N'
    if self_group in query_hits_group_iden_dict_mean:
        self_group_mean_iden = query_hits_group_iden_dict_mean[self_group]

    for_print = [query_split[1], each_extracted_split[0], '%s|%s' % (self_group, self_group_mean_iden)]
    for hit_key in query_hits_group_iden_dict_mean:
        if hit_key != self_group:
            for_print.append('%s|%s' % (hit_key, query_hits_group_iden_dict_mean[hit_key]))

    blast_hit_in_one_line_extracted_tmp1_handle.write('%s\n' % '\t'.join(for_print))
blast_hit_in_one_line_extracted_tmp1_handle.close()


# sort according to gene ID
os.system('cat %s | sort > %s' % (blast_hit_in_one_line_extracted_tmp1, blast_hit_in_one_line_extracted_tmp1_sorted))


# rename order level group ID
class_to_order_dict = {}
for each_bin in open(gtdbtk_output):
    if not each_bin.startswith('user_genome'):
        each_bin_split = each_bin.strip().split('\t')
        taxon_str = each_bin.strip().split('\t')[1]
        taxon_str_split = taxon_str.split(';')
        class_name = taxon_str_split[2]
        order_name = taxon_str_split[3]
        if order_name != 'o__':
            if class_name not in class_to_order_dict:
                class_to_order_dict[class_name] = {order_name}
            else:
                class_to_order_dict[class_name].add(order_name)

taxon_to_group_id_dict_c = {}
for taxon_c in open(taxon_to_group_id_file_c):
    taxon_c_split = taxon_c.strip().split(',')
    taxon_to_group_id_dict_c[taxon_c_split[1]] = taxon_c_split[0]

taxon_to_group_id_dict_o = {}
for taxon_o in open(taxon_to_group_id_file_o):
    taxon_o_split = taxon_o.strip().split(',')
    taxon_to_group_id_dict_o[taxon_o_split[1]] = taxon_o_split[0]

class_to_order_dict_renamed = {}
for each_class in class_to_order_dict:
    each_class_id = taxon_to_group_id_dict_c[each_class]
    orders_in_current_class = class_to_order_dict[each_class]
    orders_in_current_class_id = [taxon_to_group_id_dict_o[i] for i in orders_in_current_class]
    class_to_order_dict_renamed[each_class_id] = orders_in_current_class_id

# print(class_to_order_dict_renamed)
#
#
# all_order_list = []
# for each in class_to_order_dict_renamed:
#     for each_o in class_to_order_dict_renamed[each]:
#         all_order_list.append(each_o)
#     print('%s: %s' % (each, '\t'.join(class_to_order_dict_renamed[each])))
# print(all_order_list)

rename_dict = {'AC':'NB', 'AA':'NA', 'G':'NC', 'AH':'L0', 'Y':'M0', 'M':'J0', 'A':'H0', 'U':'KB', 'Z':'KA', 'AB':'BD', 'AG':'BC', 'F':'BA', 'W':'BB', 'AD':'BF', 'AF':'BG', 'N':'BI', 'B':'BJ', 'Q':'BH', 'C':'BE', 'L':'F0', 'T':'E0', 'S':'I0', 'D':'G0', 'I':'CE', 'X':'CC', 'H':'CG', 'R':'CH', 'V':'CF', 'AE':'CB', 'P':'CA', 'J':'CI', 'K':'CD', 'O':'O0', 'E':'D0'}


blast_hit_in_one_line_extracted_tmp1_sorted2_handle = open(blast_hit_in_one_line_extracted_tmp1_sorted2, 'w')
for each_line in open(blast_hit_in_one_line_extracted_tmp1_sorted):
    each_line_split = each_line.strip().split('\t')
    cate = each_line_split[1]
    if '_in_o' in cate:
        each_line_split_new = [each_line_split[0], each_line_split[1]]
        for iden in each_line_split[2:]:
            iden_group_old = iden.split('|')[0]
            iden_value = iden.split('|')[1]
            iden_group_new = rename_dict[iden_group_old]
            iden_renamed = '%s|%s' % (iden_group_new, iden_value)
            each_line_split_new.append(iden_renamed)
        blast_hit_in_one_line_extracted_tmp1_sorted2_handle.write('%s\n' % '\t'.join(each_line_split_new))
    else:
        blast_hit_in_one_line_extracted_tmp1_sorted2_handle.write(each_line)

blast_hit_in_one_line_extracted_tmp1_sorted2_handle.close()






