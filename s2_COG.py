import os
import numpy as np
import pandas as pd
from scipy import stats
import Kelp_HGT_config


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


########################################################## COG #########################################################

Get_COG_boxplot = '''

cd /srv/scratch/z5039045/MetaCHIP_Kelp_NM/Kelp_NM_dRep99_MetaCHIP_wd
mkdir Kelp_and_HGT_COG_func_stats
mkdir Tara_and_HGT_COG_func_stats

cp faa_files_Kelp_and_HGT_COG2014_wd/*/*_func_stats_GeneNumber.txt Kelp_and_HGT_COG_func_stats/
cp faa_files_Tara_NM_and_HGT_COG2014_wd/*/*_func_stats_GeneNumber.txt Tara_and_HGT_COG_func_stats/

# get data matrix
cd /Users/songweizhi/Desktop/Kelp_NM/COG
python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in Kelp_and_HGT_COG_func_stats -out Kelp_and_HGT_COG_func_stats_pct.txt -skip_1st_row -with_functional_description -in_percent
python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in Tara_and_HGT_COG_func_stats -out Tara_and_HGT_COG_func_stats_pct.txt -skip_1st_row -with_functional_description -in_percent

python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in Kelp_and_HGT_COG_func_stats -out Kelp_and_HGT_COG_func_stats.txt -skip_1st_row -with_functional_description
python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in Tara_and_HGT_COG_func_stats -out Tara_and_HGT_COG_func_stats.txt -skip_1st_row -with_functional_description


# get plot
cd /Users/songweizhi/Desktop/Kelp_NM/COG
Rscript ~/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i Kelp_and_HGT_COG_func_stats_pct.txt -o Kelp_and_HGT_COG_func_stats_pct.png
Rscript ~/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i Tara_and_HGT_COG_func_stats_pct.txt -o Tara_and_HGT_COG_func_stats_pct.png

'''

######################################## get details of enriched COG categories ########################################

def get_enriched_cog_id(enriched_COG_cates, HGT_annotation_results, enriched_COG_id_txt):

    enriched_COG_id_txt_handle = open(enriched_COG_id_txt, 'w')

    for enriched_COG_cate in enriched_COG_cates:
        enriched_COG_id_txt_handle.write('%s\n' % enriched_COG_cate)

        current_cate_cog_to_description_dict = {}
        current_cate_cog_to_num_dict = {}
        current_cate_identified_cog = []
        for hgt_cog in open(HGT_annotation_results):
            if not hgt_cog.startswith('Query	COG'):
                hgt_cog_split = hgt_cog.strip().split('\t')
                if len(hgt_cog_split) > 1:
                    cog_id = hgt_cog_split[1]
                    cog_cate = hgt_cog_split[2]
                    cog_description = hgt_cog_split[3]

                    if enriched_COG_cate in cog_cate:

                        # get current_cate_identified_cog
                        if cog_id not in current_cate_identified_cog:
                            current_cate_identified_cog.append(cog_id)

                        # get current_cate_cog_to_description_dict
                        current_cate_cog_to_description_dict[cog_id] = cog_description

                        # get current_cate_cog_to_num_dict
                        if cog_id not in current_cate_cog_to_num_dict:
                            current_cate_cog_to_num_dict[cog_id] = 1
                        else:
                            current_cate_cog_to_num_dict[cog_id] += 1

        for current_cate_cog in sorted(current_cate_identified_cog):
            enriched_COG_id_txt_handle.write('%s\t%s\t%s\t%s\n' % (
            enriched_COG_cate, current_cate_cog, current_cate_cog_to_num_dict[current_cate_cog],
            current_cate_cog_to_description_dict[current_cate_cog]))

    enriched_COG_id_txt_handle.close()


Kelp_HGT_annotation_results = '/Users/songweizhi/Desktop/Kelp_NM/COG/zKelp_recipient_genes_COG2014_wd/zKelp_recipient_genes_query_to_cog.txt'
Tara_HGT_annotation_results = '/Users/songweizhi/Desktop/Kelp_NM/COG/zTara_NM_recipient_genes_COG2014_wd/zTara_NM_recipient_genes_query_to_cog.txt'

Kelp_enriched_COG_cates = ['J', 'C', 'G', 'E', 'I', 'P', 'Q']
Tara_enriched_COG_cates = ['J', 'C', 'E', 'P']

Kelp_HGT_enriched_cog = '/Users/songweizhi/Desktop/Kelp_NM/COG/Kelp_and_HGT_enriched_cog.txt'
Tara_HGT_enriched_cog = '/Users/songweizhi/Desktop/Kelp_NM/COG/Tara_and_HGT_enriched_cog.txt'

get_enriched_cog_id(Kelp_enriched_COG_cates, Kelp_HGT_annotation_results, Kelp_HGT_enriched_cog)
get_enriched_cog_id(Tara_enriched_COG_cates, Tara_HGT_annotation_results, Tara_HGT_enriched_cog)


###################################### COG enrichment analysis (to be continued) #######################################

# # transpose csv file
# transpose_csv(Kelp_HGT_config.Kelp_dbCAN_df_txt, Kelp_HGT_config.Kelp_dbCAN_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Tara_dbCAN_df_txt, Kelp_HGT_config.Tara_dbCAN_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Kelp_COG_df_txt, Kelp_HGT_config.Kelp_COG_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Tara_COG_df_txt, Kelp_HGT_config.Tara_COG_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Kelp_KEGG_df_txt, Kelp_HGT_config.Kelp_KEGG_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Tara_KEGG_df_txt, Kelp_HGT_config.Tara_KEGG_df_txt_t, '\t', 0, 0)


# identified_COG_cate_list = []
# Kelp_COG_cate_list_dict = {}
# for i in open(Kelp_HGT_config.Kelp_COG_df_txt_t):
#     if not i.startswith('\t'):
#         i_split = i.strip().split('\t')
#         i_id = i_split[0]
#         i_among_MAGs_float = [float(v) for v in i_split[1:-1]]
#         absent_percent = i_among_MAGs_float.count(0)*100/len(i_among_MAGs_float)
#
#         if absent_percent <= 25:
#             Kelp_COG_cate_list_dict[i_id] = i_among_MAGs_float
#
#             if i_id not in identified_COG_cate_list:
#                 identified_COG_cate_list.append(i_id)
#
# Tara_NM_COG_cate_list_dict = {}
# for i in open(Kelp_HGT_config.Tara_COG_df_txt_t):
#     if not i.startswith('\t'):
#         i_split = i.strip().split('\t')
#         i_id = i_split[0]
#         i_among_MAGs_float = [float(v) for v in i_split[1:-1]]
#         absent_percent = i_among_MAGs_float.count(0)*100/len(i_among_MAGs_float)
#
#         if absent_percent <= 25:
#             Tara_NM_COG_cate_list_dict[i_id] = i_among_MAGs_float
#
#             if i_id not in identified_COG_cate_list:
#                 identified_COG_cate_list.append(i_id)
#
# for i in identified_COG_cate_list:
#     Kelp_MAG_i_proportion = []
#     if i in Kelp_COG_cate_list_dict:
#         Kelp_MAG_i_proportion = Kelp_COG_cate_list_dict[i]
#     Tara_MAG_i_proportion = []
#     if i in Tara_NM_COG_cate_list_dict:
#         Tara_MAG_i_proportion = Tara_NM_COG_cate_list_dict[i]
#     Kelp_MAG_i_proportion_array = np.array(Kelp_MAG_i_proportion)
#     Tara_MAG_i_proportion_array = np.array(Tara_MAG_i_proportion)

#    print(stats.normaltest(Kelp_MAG_i_proportion_array))
#    print(stats.normaltest(Tara_MAG_i_proportion_array))
    # T_test_result = stats.ttest_ind(np.array(Kelp_MAG_i_proportion), np.array(Tara_MAG_i_proportion))
    # T_test_result_pvalue = T_test_result.pvalue
    # if T_test_result_pvalue > 0.05:
    #     print(i)


# test normality in python (scipy.stats.normaltest)
# test for equal variances: scipy.stats.levene and
