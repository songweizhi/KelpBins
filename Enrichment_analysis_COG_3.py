
output_test =       '/Users/songweizhi/Desktop/COG_enrichment_analysis/COG_enrichment_test_results_Mann_Whitney_U.tab'
COG_fun_db_file =   '/Users/songweizhi/Desktop/COG_enrichment_analysis/cognames2003-2014.tab'

# file out
emriched_in_kelp =  '/Users/songweizhi/Desktop/COG_enrichment_analysis/COG_enrichment_test_emriched_in_kelp.txt'
emriched_in_tara =  '/Users/songweizhi/Desktop/COG_enrichment_analysis/COG_enrichment_test_emriched_in_tara.txt'
emriched_in_na =    '/Users/songweizhi/Desktop/COG_enrichment_analysis/COG_enrichment_test_emriched_in_na.txt'

# get cog id to function dict
cog_id_to_function_dict = {}
for each_cog in open(COG_fun_db_file):
    each_cog_split = each_cog.strip().split('\t')
    cog_id_to_function_dict[each_cog_split[0]] = each_cog_split[2]


emriched_in_kelp_handle = open(emriched_in_kelp, 'w')
emriched_in_tara_handle = open(emriched_in_tara, 'w')
emriched_in_na_handle = open(emriched_in_na, 'w')
emriched_in_kelp_handle.write('COG	P_value	Kelp	Planktonic  Mean_diff	Enriched_in Function\n')
emriched_in_tara_handle.write('COG	P_value	Kelp	Planktonic  Mean_diff	Enriched_in Function\n')
emriched_in_na_handle.write('COG	P_value	Kelp	Planktonic  Mean_diff	Enriched_in Function\n')
n = 0
for cog in open(output_test):

    if cog.startswith('COG	Kelp	Planktonic	P_value	P_value_adjusted'):
        sample_1_name = cog.strip().split('\t')[1]
        sample_2_name = cog.strip().split('\t')[2]

    else:
        cog_split = cog.strip().split('\t')
        P_value_adjusted = float(cog_split[4])

        if P_value_adjusted <= 0.05:

            mean_diff = 'NA'
            if (float(cog_split[1]) > 0) and (float(cog_split[2]) > 0):
                mean_diff = float("{0:.3f}".format(float(cog_split[1])/float(cog_split[2])))

            enriched_in = ''
            if mean_diff == 'NA':
                enriched_in = 'NA'
            elif mean_diff >= 2:
                enriched_in = sample_1_name
            elif mean_diff <= 0.5:
                enriched_in = sample_2_name

            if mean_diff == 'NA':
                emriched_in_na_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cog_split[0], P_value_adjusted, float(cog_split[1]), float(cog_split[2]), mean_diff, enriched_in, cog_id_to_function_dict[cog_split[0]]))
            elif mean_diff >= 2:
                emriched_in_kelp_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cog_split[0], P_value_adjusted, float(cog_split[1]), float(cog_split[2]), mean_diff, enriched_in, cog_id_to_function_dict[cog_split[0]]))
            elif mean_diff <= 0.5:
                emriched_in_tara_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (cog_split[0], P_value_adjusted, float(cog_split[1]), float(cog_split[2]), mean_diff, enriched_in, cog_id_to_function_dict[cog_split[0]]))

emriched_in_kelp_handle.close()
emriched_in_tara_handle.close()
emriched_in_na_handle.close()
