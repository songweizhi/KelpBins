import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from Bio import SeqIO
import matplotlib.pyplot as plt
import s0_Kelp_bins_config as KelpCfg
from s0_suspicious_HGTs import suspicious_HGTs


def get_ctg_match_cate_and_identity_distribution_plot(pwd_candidates_file_ET, pwd_iden_distribution_plot_BM, pwd_iden_distribution_plot_PG):

    # read in prediction results
    HGT_num_BM_normal = 0
    HGT_num_BM_at_end = 0
    HGT_num_BM_ctg_aln = 0
    HGT_num_PG_normal = 0
    HGT_num_PG_at_end = 0
    HGT_num_PG_ctg_aln = 0
    identity_list_BM_normal = []
    identity_list_BM_end_match = []
    identity_list_BM_full_length_match = []
    identity_list_PG_normal = []
    identity_list_PG_end_match = []
    identity_list_PG_full_length_match = []
    for each_HGT in open(pwd_candidates_file_ET):
        if not each_HGT.startswith('Gene_1'):
            each_HGT_split = each_HGT.strip().split('\t')
            identity = float(each_HGT_split[2])

            identity = 100 - identity
            end_match = each_HGT_split[4]
            full_length_match = each_HGT_split[5]
            PG_validation = each_HGT_split[6]

            # get number for normal
            if (end_match == 'no') and (full_length_match == 'no'):
                HGT_num_BM_normal += 1
                identity_list_BM_normal.append(identity)
                if PG_validation != 'NA':
                    HGT_num_PG_normal += 1
                    identity_list_PG_normal.append(identity)

            # get number for at_end
            if (end_match == 'yes') and (full_length_match == 'no'):
                HGT_num_BM_at_end += 1
                identity_list_BM_end_match.append(identity)
                if PG_validation != 'NA':
                    HGT_num_PG_at_end += 1
                    identity_list_PG_end_match.append(identity)

            # get number for Ctg_align
            if (end_match == 'no') and (full_length_match == 'yes'):
                HGT_num_BM_ctg_aln += 1
                identity_list_BM_full_length_match.append(identity)
                if PG_validation != 'NA':
                    HGT_num_PG_ctg_aln += 1
                    identity_list_PG_full_length_match.append(identity)


    #################################### plot identity distribution of identified HGTs #####################################

    ########## for BM approach ##########
    num_bins = 50
    combined_list_BM = (identity_list_BM_normal, identity_list_BM_end_match, identity_list_BM_full_length_match)
    color_list = ['g', 'orange', 'r']
    label_list = ['Normal', 'End Match', 'Full Length Match']
    plt.hist(combined_list_BM, num_bins, alpha=0.6, normed=0, stacked=1, linewidth=0, color=color_list, label=label_list, rwidth=0.85)
    lgd = plt.legend(prop={'size': 10}, ncol=1, bbox_to_anchor=(1.27, 1))

    total_HGT_num_BM = len(identity_list_BM_normal) + len(identity_list_BM_end_match) + len(identity_list_BM_full_length_match)
    #plt.title('Identity distribution of identified %s HGTs' % total_HGT_num_BM)
    plt.xlabel('Genetic variation (%)')
    plt.ylabel('Number of HGT')

    # Get plot
    plt.tight_layout()
    plt.savefig(pwd_iden_distribution_plot_BM, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()


    ########## for PG approach ##########
    num_bins = 50
    combined_list_PG = (identity_list_PG_normal, identity_list_PG_end_match, identity_list_PG_full_length_match)
    color_list = ['g', 'orange', 'r']
    label_list = ['Normal', 'End Match', 'Full Length Match']
    plt.hist(combined_list_PG, num_bins, alpha=0.6, normed=0, stacked=1, linewidth=0, color=color_list, label=label_list, rwidth=0.85)
    lgd = plt.legend(prop={'size': 10}, ncol=1, bbox_to_anchor=(1.27, 1))

    total_HGT_num_PG = len(identity_list_PG_normal) + len(identity_list_PG_end_match) + len(identity_list_PG_full_length_match)
    #plt.title('Identity distribution of identified %s HGTs' % total_HGT_num_PG)
    plt.xlabel('Genetic variation (%)')
    plt.ylabel('Number of HGT')

    # Get plot
    plt.tight_layout()
    plt.savefig(pwd_iden_distribution_plot_PG, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()


################################# plot genetic variation of PG validated HGTs (pcofg) ##################################

# get identity list
identity_list = []
gene_id_list = set()
for hgt_pcofg in open(KelpCfg.HGT_PG_validated_txt_pcofg):
    if not hgt_pcofg.startswith('Gene_1'):
        hgt_pcofg_split = hgt_pcofg.strip().split('\t')
        concatenated_genes = '%s___%s' % (hgt_pcofg_split[0], hgt_pcofg_split[1])
        identity = float(hgt_pcofg_split[2])

        # only plot those with high confidence
        if concatenated_genes not in suspicious_HGTs:
            identity_list.append(identity)
            gene_id_list.add(hgt_pcofg_split[0])
            gene_id_list.add(hgt_pcofg_split[1])


# transfer identity to genetic variation
genetic_variation_list = []
for each in identity_list:
    genetic_variation_list.append(100 -each)


# Get plot
num_bins = 35
plt.hist(genetic_variation_list, num_bins, alpha=0.6, normed=0, linewidth=0, color='g', rwidth=0.8)
#plt.title('Genetic variation of identified HGTs')
plt.xlabel('Genetic variation (%)')
plt.ylabel('Number of HGT')
plt.tight_layout()
plt.savefig(KelpCfg.pwd_plot_genetic_variation_PG_validated_pcofg, bbox_inches='tight', dpi=300)
plt.close()
plt.clf()


get_ctg_match_cate_and_identity_distribution_plot(KelpCfg.HGT_PG_txt_pcofg, KelpCfg.pwd_plot_genetic_variation_BM_pcofg, KelpCfg.pwd_plot_genetic_variation_PG_pcofg)


################################# get the number of HGT at different variations levels #################################

range_id = ['0-5', '5-10', '10-15', '15-20', '20-25', '25-30', '30-35']

genetic_variation_group_dict = {}
num_0 = 0
num_5 = 0
num_10 = 0
num_15 = 0
num_20 = 0
num_25 = 0
num_30 = 0
for each_iden in identity_list:
    genetic_variation = float("{0:.2f}".format(100 - each_iden))

    if genetic_variation <= 5:
        if '0-5' not in genetic_variation_group_dict:
            genetic_variation_group_dict['0-5'] = [genetic_variation]
        else:
            genetic_variation_group_dict['0-5'].append(genetic_variation)

    if 5 < genetic_variation <= 10:
        if '5-10' not in genetic_variation_group_dict:
            genetic_variation_group_dict['5-10'] = [genetic_variation]
        else:
            genetic_variation_group_dict['5-10'].append(genetic_variation)

    if 10 < genetic_variation <= 15:
        if '10-15' not in genetic_variation_group_dict:
            genetic_variation_group_dict['10-15'] = [genetic_variation]
        else:
            genetic_variation_group_dict['10-15'].append(genetic_variation)

    if 15 < genetic_variation <= 20:
        if '15-20' not in genetic_variation_group_dict:
            genetic_variation_group_dict['15-20'] = [genetic_variation]
        else:
            genetic_variation_group_dict['15-20'].append(genetic_variation)

    if 20 < genetic_variation <= 25:
        if '20-25' not in genetic_variation_group_dict:
            genetic_variation_group_dict['20-25'] = [genetic_variation]
        else:
            genetic_variation_group_dict['20-25'].append(genetic_variation)

    if 25 < genetic_variation <= 30:
        if '25-30' not in genetic_variation_group_dict:
            genetic_variation_group_dict['25-30'] = [genetic_variation]
        else:
            genetic_variation_group_dict['25-30'].append(genetic_variation)

    if 30 < genetic_variation <= 35:
        if '30-35' not in genetic_variation_group_dict:
            genetic_variation_group_dict['30-35'] = [genetic_variation]
        else:
            genetic_variation_group_dict['30-35'].append(genetic_variation)

range_num_handle = open(KelpCfg.pwd_genetic_variation_PG_validated_pcofg_range_num, 'w')
range_num_handle.write('Range\tNumber\tPercent\n')
for each_range in range_id:
    each_range_num = len(genetic_variation_group_dict[each_range])
    each_range_percent = float("{0:.2f}".format(each_range_num*100/len(identity_list)))
    range_num_handle.write('%s\t%s\t%s\n' % (each_range, each_range_num, each_range_percent))
range_num_handle.close()



#####

num_0 = 0
num_5 = 0
num_10 = 0
num_15 = 0
num_20 = 0
num_25 = 0
num_30 = 0
for each_iden in identity_list:
    gv = 100 - each_iden

    if gv <= 2.5:
        num_0 += 1
    elif 2.5 < gv <= 7.5:
        num_5 += 1
    elif 7.5 < gv <= 12.5:
        num_10 += 1
    elif 12.5 < gv <= 17.5:
        num_15 += 1
    elif 17.5 < gv <= 22.5:
        num_20 += 1
    elif 22.5 < gv <= 27.5:
        num_25 += 1
    elif 27.5 < gv <= 32.5:
        num_30 += 1


total_num = num_0 + num_5 + num_10 + num_15 + num_20 + num_25 + num_30
print(total_num)
print('0\t5\t10\t15\t20\t25\t30')
print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (num_0, num_5, num_10, num_15, num_20, num_25, num_30))
print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (float("{0:.2f}".format(num_0*100/total_num)),
                                      float("{0:.2f}".format(num_5*100/total_num)),
                                      float("{0:.2f}".format(num_10*100/total_num)),
                                      float("{0:.2f}".format(num_15*100/total_num)),
                                      float("{0:.2f}".format(num_20*100/total_num)),
                                      float("{0:.2f}".format(num_25*100/total_num)),
                                      float("{0:.2f}".format(num_30*100/total_num))))


################################ extract sequencei of normal HGTs after manually check #################################

pwd_HGT_pcofg_faa_manually_checked_file_handle = open(KelpCfg.pwd_HGT_pcofg_faa_manually_checked, 'w')
for each_seq in SeqIO.parse(KelpCfg.all_combined_faa, 'fasta'):
    if each_seq.id in gene_id_list:
        pwd_HGT_pcofg_faa_manually_checked_file_handle.write('>%s\n' % each_seq.id)
        pwd_HGT_pcofg_faa_manually_checked_file_handle.write('%s\n' % str(each_seq.seq))
pwd_HGT_pcofg_faa_manually_checked_file_handle.close()


# ################################ plot number of HGT detected from each genome and group ################################
#
# # store genome size in dict
# genome_size_dict = {}
# for each_size in open(KelpCfg.genome_size_file):
#     if each_size.strip() != 'Genome\tSize(Mbp)':
#         each_size_split = each_size.strip().split('\t')
#         genome_file_name = each_size_split[0]
#         genome_file_name_no_ext = '.'.join(genome_file_name.split('.')[:-1])
#         genome_size_Mbp = float(each_size_split[1])
#         genome_size_dict[genome_file_name_no_ext] = genome_size_Mbp
#
# # get genome_to_HgtNum_dict
# genome_to_HgtNum_dict = {}
# for each_pair in open(KelpCfg.HGT_PG_validated_txt_pcofg):
#     if not each_pair.startswith('Gene_1'):
#
#         each_pair_split = each_pair.strip().split('\t')
#         gene_1 = each_pair_split[0]
#         gene_2 = each_pair_split[1]
#         genome_1 = '_'.join(gene_1.split('_')[:-1])
#         genome_2 = '_'.join(gene_2.split('_')[:-1])
#
#         concatenated_genes = '%s___%s' % (gene_1, gene_2)
#
#         # only plot those with high confidence
#         if concatenated_genes not in suspicious_HGTs:
#
#             if genome_1 not in genome_to_HgtNum_dict:
#                 genome_to_HgtNum_dict[genome_1] = 1
#             else:
#                 genome_to_HgtNum_dict[genome_1] += 1
#
#             if genome_2 not in genome_to_HgtNum_dict:
#                 genome_to_HgtNum_dict[genome_2] = 1
#             else:
#                 genome_to_HgtNum_dict[genome_2] += 1
#
#
#
#
#
#
# # define color list
# color_list = ['black', 'silver', 'darksalmon', 'blueviolet', 'sienna', 'tan', 'purple', 'gold', 'palegreen',
#               'paleturquoise', 'slategray', 'royalblue', 'plum', 'olivedrab', 'seagreen', 'darkorchid',
#               'darkkhaki'] * 100
#
# output_txt_handle = open(pwd_candidates_file_ET_validated_STAT_genome_txt, 'w')
# output_txt_handle.write('Group\tGenome\tSize\tHGT\tHGT/Mbp\n')
# genome_list_according_group = []
# genome_list_according_group_with_group = []
# num_list_according_group = []
# num_list_according_group_norm = []
# color_list_according_group = []
# color_index = 1
# for each_group_id in sorted(group_id_list):
#     current_group_genomes = group_to_genome_dict[each_group_id]
#     for each_genome in current_group_genomes:
#         each_genome_with_group = '(%s)%s' % (each_group_id, each_genome)
#         if each_genome in genome_to_HgtNum_dict:
#             genome_list_according_group.append(each_genome)
#             genome_list_according_group_with_group.append(each_genome_with_group)
#
#             num_HGT = genome_to_HgtNum_dict[each_genome]
#             num_HGT_norm = num_HGT / genome_size_dict[each_genome]
#             num_HGT_norm = float("{0:.2f}".format(num_HGT_norm))
#
#             num_list_according_group.append(num_HGT)
#             num_list_according_group_norm.append(num_HGT_norm)
#
#             color_list_according_group.append(color_list[color_index])
#             for_out = '%s\t%s\t%s\t%s\t%s\n' % (
#             each_group_id, each_genome, genome_size_dict[each_genome], num_HGT, num_HGT_norm)
#             output_txt_handle.write(for_out)
#     color_index += 1
# output_txt_handle.close()
#
# # write stats to file
# HGT_PG_STAT_handle = open(pwd_candidates_file_ET_validated_STAT_group_txt, 'w')
# n = 0
#
# if grouping_level == 'x':
#     HGT_PG_STAT_handle.write('Group\tSize(Mbp)\tHGT\tHGT/Mbp\n')
# else:
#     HGT_PG_STAT_handle.write('Group\tSize(Mbp)\tHGT\tHGT/Mbp\tTaxon\n')
#
# for each_g in group_list_uniq:
#
#     if grouping_level == 'x':
#         for_out = '%s\t%s\t%s\t%s\n' % (each_g, float("{0:.3f}".format(group_to_length_dict[each_g])), group_list_uniq_count[n],group_list_uniq_count_normalized[n])
#     else:
#         for_out = '%s\t%s\t%s\t%s\t%s\n' % (each_g, float("{0:.3f}".format(group_to_length_dict[each_g])), group_list_uniq_count[n],group_list_uniq_count_normalized[n], group_2_taxon_dict[each_g])
#
#     HGT_PG_STAT_handle.write(for_out)
#     n += 1
# HGT_PG_STAT_handle.close()
#
#
# # set xticks fontsize for genome plot
# xticks_fontsize_genome = 8
# if 25 < len(genome_list_according_group_with_group) <= 50:
#     xticks_fontsize_genome = 5
# elif 50 < len(genome_list_according_group_with_group) <= 100:
#     xticks_fontsize_genome = 4
# elif 100 < len(genome_list_according_group_with_group) <= 500:
#     xticks_fontsize_genome = 3
# elif len(genome_list_according_group_with_group) > 500:
#     xticks_fontsize_genome = 2
#
#
# # set figure size
# plt.figure(figsize=(20, 20))
# x_range_genome = range(len(num_list_according_group))
#
# # subplot 1
# plt.subplot(121)
# plt.bar(x_range_genome, num_list_according_group, alpha=0.5, linewidth=0, color=color_list_according_group)
# plt.xticks(x_range_genome, genome_list_according_group_with_group, rotation=315, fontsize=xticks_fontsize_genome,
#            horizontalalignment='left')
# plt.ylabel('Number of HGT')
#
#
# # subplot 3
# plt.subplot(122)
# plt.bar(x_range_genome, num_list_according_group_norm, alpha=0.5, linewidth=0, color=color_list_according_group)
# plt.xticks(x_range_genome, genome_list_according_group_with_group, rotation=315, fontsize=xticks_fontsize_genome,
#            horizontalalignment='left')
# plt.xlabel('Genome')
# plt.ylabel('Number of HGT / Mbp sequences')
#
#
# # plot layout
# plt.subplots_adjust(wspace=0.2, top=0.95)
#
# # save plot
# plt.tight_layout()
# plt.savefig(pwd_candidates_file_ET_validated_STAT_png, dpi=300)
#
#
#
#
#
#

