import os
import glob
import scipy
from scipy import stats
import pandas as pd
from pandas import DataFrame
from matplotlib import pyplot

cogs_enriched_in_kelp = '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/COG_enrichment_test_emriched_in_kelp.txt'
query_to_cog_file_folder = '/Users/songweizhi/Desktop/Kelp_NM/Figure_data/Figure_2_COG_stats_absolute/faa_files_Kelp_query_to_cog_files'
query_to_cog_file_re = '%s/*_query_to_cog.txt' % query_to_cog_file_folder
COG_fun_file =              '/Users/songweizhi/DB/COG2014/fun2003-2014.tab'

COG_order =                 'JAKLBDYVTMNZWUOCGEFHIPQRS'


COG_cate_to_fun_dict = {}
for COG_cate in open(COG_fun_file):
    COG_cate_split = COG_cate.strip().split('\t')
    COG_cate_to_fun_dict[COG_cate_split[0]] = COG_cate_split[1]


query_to_cog_file_list = [os.path.basename(file_name) for file_name in glob.glob(query_to_cog_file_re)]


cate_to_identified_cog_dict = {}
for query_to_cog_file in query_to_cog_file_list:
    pwd_query_to_cog_file = '%s/%s' % (query_to_cog_file_folder, query_to_cog_file)
    for query_to_cog in open(pwd_query_to_cog_file):
        if not query_to_cog.startswith('Query	COG	Category	Description'):
            query_to_cog_split = query_to_cog.strip().split('\t')
            if len(query_to_cog_split) > 1:
                cog_id = query_to_cog_split[1]
                cog_cate = query_to_cog_split[2]
                for each_cate in cog_cate:
                    if each_cate not in cate_to_identified_cog_dict:
                        cate_to_identified_cog_dict[each_cate] = [cog_id]
                    else:
                        if cog_id not in cate_to_identified_cog_dict[each_cate]:
                            cate_to_identified_cog_dict[each_cate].append(cog_id)


kelp_enriched_cate_to_cog_dict = {}
for kelp_enriched_cog in open(cogs_enriched_in_kelp):
    if not kelp_enriched_cog.startswith('Category	COG	P_value	Kelp'):
        kelp_enriched_cog_split = kelp_enriched_cog.strip().split('\t')
        kelp_enriched_cog_cate = kelp_enriched_cog_split[0]
        kelp_enriched_cog_id = kelp_enriched_cog_split[1]
        for each_kelp_enriched_cog_cate in kelp_enriched_cog_cate:
            if each_kelp_enriched_cog_cate not in kelp_enriched_cate_to_cog_dict:
                kelp_enriched_cate_to_cog_dict[each_kelp_enriched_cog_cate] = [kelp_enriched_cog_id]
            else:
                if kelp_enriched_cog_id not in kelp_enriched_cate_to_cog_dict[each_kelp_enriched_cog_cate]:
                    kelp_enriched_cate_to_cog_dict[each_kelp_enriched_cog_cate].append(kelp_enriched_cog_id)


print('COG_cate\tAll\tEnriched')
for COG_cate in COG_order:

    COG_cate_id_num = 0
    if COG_cate in cate_to_identified_cog_dict:
        COG_cate_id_num = len(cate_to_identified_cog_dict[COG_cate])

    COG_cate_id_num_enriched = 0
    if COG_cate in kelp_enriched_cate_to_cog_dict:
        COG_cate_id_num_enriched = len(kelp_enriched_cate_to_cog_dict[COG_cate])

    print('%s_%s\t%s\t%s' % (COG_cate, COG_cate_to_fun_dict[COG_cate], COG_cate_id_num, COG_cate_id_num_enriched))
