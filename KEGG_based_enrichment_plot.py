import os
import numpy as np
import pandas as pd
from scipy import stats


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


# file in
ko_level    = 'B'
mag_source  = 'Tara_NM'  # Kelp
min_HGT_num = 10

faa_files_Kelp_and_HGT_B_GeneNumber                 = '/Users/songweizhi/Desktop/Kelp_NM/KEGG/faa_files_%s_and_HGT_%s_GeneNumber.txt'               % (mag_source, ko_level)
faa_files_Kelp_and_HGT_B_GeneNumber_pct             = '/Users/songweizhi/Desktop/Kelp_NM/KEGG/faa_files_%s_and_HGT_%s_GeneNumber_pct.txt'           % (mag_source, ko_level)
KEGG_DB_ko                                          = '/Users/songweizhi/DB/KEGG_2016-09-26/ko00001.keg'

# file out
faa_files_Kelp_and_HGT_B_GeneNumber_t               = '/Users/songweizhi/Desktop/Kelp_NM/KEGG/faa_files_%s_and_HGT_%s_GeneNumber_t.txt'             % (mag_source, ko_level)
faa_files_Kelp_and_HGT_B_GeneNumber_pct_t           = '/Users/songweizhi/Desktop/Kelp_NM/KEGG/faa_files_%s_and_HGT_%s_GeneNumber_pct_t.txt'         % (mag_source, ko_level)
faa_files_Kelp_and_HGT_B_GeneNumber_pct_t_minHGT    = '/Users/songweizhi/Desktop/Kelp_NM/KEGG/faa_files_%s_and_HGT_%s_GeneNumber_pct_t_min%s.txt'   % (mag_source, ko_level, min_HGT_num)
faa_files_Kelp_and_HGT_B_GeneNumber_pct_minHGT      = '/Users/songweizhi/Desktop/Kelp_NM/KEGG/faa_files_%s_and_HGT_%s_GeneNumber_pct_min%s.txt'     % (mag_source, ko_level, min_HGT_num)


# store ko functions in dict
As_description_dict = {}
Bs_description_dict = {}
Cs_description_dict = {}
Ds_description_dict = {}
D2ABCD_dict = {}
current_A = ''
current_B = ''
current_C = ''
for each_line in open(KEGG_DB_ko):
    if each_line[0] in ['A', 'B', 'C', 'D']:
        each_line_split = each_line.strip().split(' ')

        if each_line[0] == 'A':
            current_A_id = each_line_split[0]
            current_A_description = ' '.join(each_line_split[1:])
            current_A = current_A_id
            As_description_dict[current_A_id] = current_A_description

        elif each_line[0] == 'B':
            if len(each_line_split) > 1:
                current_B_id = each_line_split[2]
                current_B_description = ' '.join(each_line_split[3:])
                current_B = current_B_id
                Bs_description_dict[current_B_id] = current_B_description

        elif each_line[0] == 'C':
            current_C_id = each_line_split[4]
            current_C_description = ' '.join(each_line_split[5:])
            current_C = current_C_id
            Cs_description_dict[current_C_id] = current_C_description

        elif each_line[0] == 'D':
            current_D_id = each_line_split[6]
            current_D_description = ' '.join(each_line_split[7:])
            Ds_description_dict[current_D_id] = current_D_description
            ABCD_value = 'A_%s|B_%s|C_%s|D_%s' % (current_A, current_B, current_C, current_D_id)
            if current_D_id not in D2ABCD_dict:
                D2ABCD_dict[current_D_id] = [ABCD_value]
            elif (current_D_id in D2ABCD_dict) and (ABCD_value not in D2ABCD_dict[current_D_id]):
                D2ABCD_dict[current_D_id].append(ABCD_value)

# transpose csv file
transpose_csv(faa_files_Kelp_and_HGT_B_GeneNumber, faa_files_Kelp_and_HGT_B_GeneNumber_t, '\t', 0, 0)
transpose_csv(faa_files_Kelp_and_HGT_B_GeneNumber_pct, faa_files_Kelp_and_HGT_B_GeneNumber_pct_t, '\t', 0, 0)


# get ko_to_hgt_num_dict
ko_to_hgt_num_dict = {}
for each_line in open(faa_files_Kelp_and_HGT_B_GeneNumber_t):
    if not each_line.startswith('\t'):
        each_line_split = each_line.strip().split('\t')
        function_id = each_line_split[0]
        num_in_hgt = int(each_line_split[-1])
        ko_to_hgt_num_dict[function_id] = num_in_hgt

os.system('rm %s' % faa_files_Kelp_and_HGT_B_GeneNumber_t)


faa_files_Kelp_and_HGT_B_GeneNumber_pct_t_minHGT_handle = open(faa_files_Kelp_and_HGT_B_GeneNumber_pct_t_minHGT, 'w')
for each_line_pct in open(faa_files_Kelp_and_HGT_B_GeneNumber_pct_t):
    if each_line_pct.startswith('\t'):
        faa_files_Kelp_and_HGT_B_GeneNumber_pct_t_minHGT_handle.write(each_line_pct)
    else:
        each_line_pct_split = each_line_pct.strip().split('\t')
        function_id = each_line_pct_split[0]
        pct_among_mags = each_line_pct_split[1:-1]
        pct_in_hgt = each_line_pct_split[-1]

        function_id_absolute_num = ko_to_hgt_num_dict[function_id]

        function_id_with_description = ''
        if ko_level == 'B':
            function_id_with_description = '%s_%s' % (function_id, Bs_description_dict[function_id])
        if ko_level == 'C':
            function_id_with_description = '%s_%s' % (function_id, Cs_description_dict[function_id])

        if function_id_absolute_num >= min_HGT_num:
            faa_files_Kelp_and_HGT_B_GeneNumber_pct_t_minHGT_handle.write('%s\t%s\t%s\n' % (function_id_with_description, '\t'.join(pct_among_mags), pct_in_hgt))

faa_files_Kelp_and_HGT_B_GeneNumber_pct_t_minHGT_handle.close()


# transpose back
transpose_csv(faa_files_Kelp_and_HGT_B_GeneNumber_pct_t_minHGT, faa_files_Kelp_and_HGT_B_GeneNumber_pct_minHGT, '\t', 0, 0)

os.system('rm %s' % faa_files_Kelp_and_HGT_B_GeneNumber_pct_t)
os.system('rm %s' % faa_files_Kelp_and_HGT_B_GeneNumber_pct_t_minHGT)



Get_KEGG_boxplot = '''

cd /Users/songweizhi/Desktop/Kelp_NM/KEGG

Rscript ~/PycharmProjects/BioSAK/BioSAK/KEGG_boxplot_last1row.R -i faa_files_Kelp_and_HGT_B_GeneNumber_pct_min10.txt -o faa_files_Kelp_and_HGT_B_GeneNumber_pct_min10.png
Rscript ~/PycharmProjects/BioSAK/BioSAK/KEGG_boxplot_last1row.R -i faa_files_Kelp_and_HGT_C_GeneNumber_pct_min10.txt -o faa_files_Kelp_and_HGT_C_GeneNumber_pct_min10.png

Rscript ~/PycharmProjects/BioSAK/BioSAK/KEGG_boxplot_last1row.R -i faa_files_Tara_NM_and_HGT_B_GeneNumber_pct_min10.txt -o faa_files_Tara_NM_and_HGT_B_GeneNumber_pct_min10.png
Rscript ~/PycharmProjects/BioSAK/BioSAK/KEGG_boxplot_last1row.R -i faa_files_Tara_NM_and_HGT_C_GeneNumber_pct_min10.txt -o faa_files_Tara_NM_and_HGT_C_GeneNumber_pct_min10.png

'''
