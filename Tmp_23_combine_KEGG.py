import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)
    return folder_list


def get_bin_to_ko_percent_dict(pwd_ko_stats, ko_id_all):

    bin_to_ko_percent = {}
    for each_ko in open(pwd_ko_stats):
        if not each_ko.startswith('Level'):
            each_ko_split = each_ko.strip().split('\t')
            each_ko_id = each_ko_split[1]
            each_ko_percent = float(each_ko_split[3])
            if each_ko_id not in ko_id_all:
                ko_id_all.append(each_ko_id)

            bin_to_ko_percent[each_ko_id] = each_ko_percent

    return bin_to_ko_percent


KEGG_annot_results_folder = '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_all_KEGG'
KEGG_wd_folder_list = get_no_hidden_folder_list(KEGG_annot_results_folder)

bin_id_list = []
ko_id_B_all = []
ko_id_C_all = []
ko_percent_all_bins_B = {}
ko_percent_all_bins_C = {}
for each_folder in KEGG_wd_folder_list:

    bin_id = each_folder.split('_KEGG_wd')[0]
    bin_id_list.append(bin_id)

    pwd_ko_stats_A = '%s/%s/%s_summary_level_A.txt' % (KEGG_annot_results_folder, each_folder, bin_id)
    pwd_ko_stats_B = '%s/%s/%s_summary_level_B.txt' % (KEGG_annot_results_folder, each_folder, bin_id)
    pwd_ko_stats_C = '%s/%s/%s_summary_level_C.txt' % (KEGG_annot_results_folder, each_folder, bin_id)
    pwd_ko_stats_D = '%s/%s/%s_summary_level_D.txt' % (KEGG_annot_results_folder, each_folder, bin_id)

    bin_to_ko_percent_B = get_bin_to_ko_percent_dict(pwd_ko_stats_B, ko_id_B_all)
    bin_to_ko_percent_C = get_bin_to_ko_percent_dict(pwd_ko_stats_C, ko_id_C_all)

    ko_percent_all_bins_B[bin_id] = bin_to_ko_percent_B
    ko_percent_all_bins_C[bin_id] = bin_to_ko_percent_C


bin_id_list = sorted(bin_id_list)
ko_id_B_all = sorted(ko_id_B_all)
ko_id_C_all = sorted(ko_id_C_all)


ko_percent_lol_B = []
ko_percent_lol_C = []
for each_bin in bin_id_list:
    current_bin_ko_B_percent_dict = ko_percent_all_bins_B[each_bin]
    current_bin_ko_C_percent_dict = ko_percent_all_bins_C[each_bin]
    percent_list_B = []
    for ko_id_B in ko_id_B_all:
        ko_B_percent = 0
        if ko_id_B in current_bin_ko_B_percent_dict:
            ko_B_percent = current_bin_ko_B_percent_dict[ko_id_B]
        percent_list_B.append(ko_B_percent)

    percent_list_C = []
    for ko_id_C in ko_id_C_all:
        ko_C_percent = 0
        if ko_id_C in current_bin_ko_C_percent_dict:
            ko_C_percent = current_bin_ko_C_percent_dict[ko_id_C]
        percent_list_C.append(ko_C_percent)

    ko_percent_lol_B.append(percent_list_B)
    ko_percent_lol_C.append(percent_list_C)

    # percent_list_str_B = [str(i) for i in percent_list_B]
    # percent_list_str_C = [str(i) for i in percent_list_C]
    # for_out_B = '%s,%s' % (each_bin, ','.join(percent_list_str_B))
    # for_out_C = '%s,%s' % (each_bin, ','.join(percent_list_str_C))
    # print(for_out_B)
    # print(for_out_C)


ko_B_percent_df = pd.DataFrame(np.array(ko_percent_lol_B), index=bin_id_list, columns=ko_id_B_all)
ko_C_percent_df = pd.DataFrame(np.array(ko_percent_lol_C), index=bin_id_list, columns=ko_id_C_all)



print(ko_B_percent_df.shape)


# print(ko_B_percent_df.shape)
#
# print(ko_B_percent_df.index)
#
#
# for each in ko_B_percent_df.index:
#     print(each)
#

#print(ko_B_percent_df)



# ko_B_percent_df_no_NA = ko_B_percent_df.loc[:, ko_B_percent_df.columns != 'NA']
# ko_C_percent_df_no_NA = ko_C_percent_df.loc[:, ko_C_percent_df.columns != 'NA']
#
# ko_B_percent_df_NA = ko_B_percent_df.loc[:, ko_B_percent_df.columns == 'NA']
#
# print(ko_B_percent_df_NA)

# print(ko_percent_df.index)
# print(ko_percent_df.columns)
# print(ko_percent_df.values)
# print(ko_C_percent_df.shape)
# print(ko_C_percent_df)


# ko_C_percent_df_no_NA.boxplot(grid=False, figsize=(100, 30))
# plt.savefig('/Users/songweizhi/Desktop/test_C.png', dpi=300)


# ko_B_percent_df_no_NA.boxplot(grid=False, figsize=(25, 15))


# HGT_value = ko_B_percent_df.loc['HGT_pcofg'].tolist()

# print(HGT_value)
# print(ko_id_B_all)
# print(len(HGT_value))
# print(len(ko_id_B_all))

# marker = [u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o']

# color_list = ['r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r']

# color_list = 'rrrrrrrrrrrrrrrrbbbrrrrrbbbrrrrrrrrrrggggrrrrrrrrrrr'
# print(len(color_list))

# print(HGT_value)
# plt.plot(HGT_value, linewidth=0, marker='o', markersize=7, markeredgewidth=0, color=color_list)
# plt.plot(HGT_value, 'ro')


# #plt.scatter(HGT_value, ko_id_B_all)
# plt.xticks(rotation=270)
# plt.savefig('/Users/songweizhi/Desktop/test_B.png', dpi=300)


# #print(len(ko_B_percent_df_NA.values))

