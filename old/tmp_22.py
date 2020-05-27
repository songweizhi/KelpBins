import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def get_ctg_match_cate_and_identity_distribution_plot(pwd_candidates_file_ET, pwd_plot_ctg_match_cate, pwd_iden_distribution_plot_BM, pwd_iden_distribution_plot_PG, pwd_plot_identity_distribution_PG_normal):

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


    ################################################## plot at_ends_stat ###################################################

    n_groups = 2
    normal_list = (HGT_num_BM_normal, HGT_num_PG_normal)
    at_end_list = (HGT_num_BM_at_end, HGT_num_PG_at_end)
    ctg_match_list = (HGT_num_BM_ctg_aln, HGT_num_PG_ctg_aln)

    # create plot
    index = np.arange(n_groups)
    bar_width = 0.15

    plt.bar(index, normal_list, bar_width, alpha=0.4, color='g', label='Normal', align='center')
    plt.bar(index + bar_width, at_end_list, bar_width, alpha=0.4, color='orange', label='End Match', align='center')
    plt.bar(index + bar_width * 2, ctg_match_list, bar_width, alpha=0.4, color='r', label='Full Length Match', align='center')

    plt.ylabel('Number of predicted HGTs')
    plt.title('Location of predicted HGTs')
    plt.xticks(index + 1.5 * bar_width, ('Best-match', 'Phylogenetic'))
    lgd = plt.legend(prop={'size': 10}, ncol=1, bbox_to_anchor=(1.27, 1))

    plt.tight_layout()
    plt.savefig(pwd_plot_ctg_match_cate, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.close()


    #################################### plot identity distribution of identified HGTs #####################################

    ########## for BM approach ##########
    num_bins = 50
    combined_list_BM = (identity_list_BM_normal, identity_list_BM_end_match, identity_list_BM_full_length_match)
    color_list = ['g', 'orange', 'r']
    label_list = ['Normal', 'End Match', 'Full Length Match']
    plt.hist(combined_list_BM, num_bins, alpha=0.6, normed=0, stacked=1, linewidth=0, color=color_list, label=label_list, rwidth=0.85)
    lgd = plt.legend(prop={'size': 10}, ncol=1, bbox_to_anchor=(1.27, 1))

    total_HGT_num_BM = len(identity_list_BM_normal) + len(identity_list_BM_end_match) + len(identity_list_BM_full_length_match)
    plt.title('Identity distribution of identified %s HGTs' % total_HGT_num_BM)
    plt.xlabel('Identity (%)')
    plt.ylabel('Number of identified HGT')

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
    plt.title('Identity distribution of identified %s HGTs' % total_HGT_num_PG)
    plt.xlabel('Identity (%)')
    plt.ylabel('Number of identified HGT')

    # Get plot
    plt.tight_layout()
    plt.savefig(pwd_iden_distribution_plot_PG, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()


    ########## for PG approach (normal) ##########
    num_bins = 50
    plt.hist(identity_list_PG_normal, num_bins, alpha=0.6, normed=0, linewidth=0, color='g', rwidth=0.85)

    plt.title('Identity distribution of identified HGTs')
    plt.xlabel('Identity (%)')
    plt.ylabel('Number of identified HGT')

    # Get plot
    plt.tight_layout()
    plt.savefig(pwd_plot_identity_distribution_PG_normal, bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()


wd = '/Users/songweizhi/Desktop/555555/TT_90MGs'
os.chdir(wd)

pwd_candidates_file_ET = 'TT_90MGs_PG.txt'
pwd_plot_at_ends_number = 'TT_90MGs_plot_ctg_match_category'
pwd_plot_identity_distribution_BM = 'TT_90MGs_plot_HGT_identity_BM.png'
pwd_plot_identity_distribution_PG = 'TT_90MGs_p10_plot_HGT_identity_PG.png'
pwd_plot_identity_distribution_PG_normal = 'TT_90MGs_p10_plot_HGT_identity_PG_normal.png'

#get_ctg_match_cate_and_identity_distribution_plot(pwd_candidates_file_ET, pwd_plot_at_ends_number, pwd_plot_identity_distribution_BM, pwd_plot_identity_distribution_PG, pwd_plot_identity_distribution_PG_normal)

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
        end_match = each_HGT_split[4]
        full_length_match = each_HGT_split[5]
        PG_validation = each_HGT_split[6]
