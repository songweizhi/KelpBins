import os
import numpy as np
import pandas as pd


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


def get_suppl_table(ko_num_df_t, ko_pct_df_t, KO_description_dict, ko_suppl_txt, ko_level):

    ko_to_HGT_num_dict = {}
    for ko_num in open(ko_num_df_t):
        if not ko_num.startswith('\t'):
            ko_num_split = ko_num.strip().split('\t')
            ko_id = ko_num_split[0]
            ko_hgt_num = int(ko_num_split[-1])
            ko_to_HGT_num_dict[ko_id] = ko_hgt_num

    ko_suppl_txt_handle = open(ko_suppl_txt, 'w')
    ko_suppl_txt_handle.write('KO_%s\tMAG_mean(%s)\tHGT(%s)\tHGT_num\tPct_diff\tDescription\n' % (ko_level, '%', '%'))
    for ko_pct in open(ko_pct_df_t):
        if not ko_pct.startswith('\t'):
            ko_pct_split = ko_pct.strip().split('\t')
            ko_id = ko_pct_split[0]
            ko_pct_MAG = ko_pct_split[1:-1]
            ko_pct_MAG_float = [float(i) for i in ko_pct_MAG]
            ko_pct_HGT = float(ko_pct_split[-1])
            ko_pct_MAG_float_mean = float("{0:.2f}".format(np.mean(np.array(ko_pct_MAG_float))))
            ko_hgt_num = ko_to_HGT_num_dict[ko_id]

            if ko_hgt_num != 0:
                ko_description = KO_description_dict[ko_id]
                pct_diff = float("{0:.2f}".format(ko_pct_HGT/np.mean(np.array(ko_pct_MAG_float))))
                ko_suppl_txt_handle.write('%s_%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_level, ko_id, ko_pct_MAG_float_mean, ko_pct_HGT, ko_hgt_num, pct_diff, ko_description))

    ko_suppl_txt_handle.close()


def get_ko2description_dict(ko00001_keg):

    ko_list_C = []
    As_description_dict = {}
    Bs_description_dict = {}
    Cs_description_dict = {}
    Ds_description_dict = {}
    D2ABCD_dict = {}

    current_A = ''
    current_B = ''
    current_C = ''
    for each_line in open(ko00001_keg):
        if each_line[0] in ['A', 'B', 'C', 'D']:
            each_line_split = each_line.strip().split(' ')

            if each_line[0] == 'A':
                current_A_id = each_line_split[0]
                current_A_description = ' '.join(each_line_split[1:-1])
                current_A = current_A_id
                As_description_dict[current_A_id] = current_A_description

            elif each_line[0] == 'B':
                if len(each_line_split) > 1:

                    print(each_line_split)
                    current_B_id = each_line_split[2]
                    current_B_description = ' '.join(each_line_split[3:])
                    current_B = current_B_id
                    Bs_description_dict[current_B_id] = current_B_description

            elif each_line[0] == 'C':
                current_C_id = each_line_split[4]
                current_C_description = ' '.join(each_line_split[5:])
                current_C = current_C_id
                Cs_description_dict[current_C_id] = current_C_description
                ko_list_C.append(current_C_id)

            elif each_line[0] == 'D':
                current_D_id = each_line_split[6]
                current_D_description = ' '.join(each_line_split[7:-1])
                Ds_description_dict[current_D_id] = current_D_description
                ABCD_value = 'A_%s|B_%s|C_%s|D_%s' % (current_A, current_B, current_C, current_D_id)
                if current_D_id not in D2ABCD_dict:
                    D2ABCD_dict[current_D_id] = [ABCD_value]
                elif (current_D_id in D2ABCD_dict) and (ABCD_value not in D2ABCD_dict[current_D_id]):
                    D2ABCD_dict[current_D_id].append(ABCD_value)

    return As_description_dict, Bs_description_dict, Cs_description_dict, Ds_description_dict, D2ABCD_dict, ko_list_C


############################################## file I/O (dRep95) #######################################################

wd = '/Users/songweizhi/Desktop/Kelp_NM/KEGG'
os.chdir(wd)

# input files
Kelp_ko_B_num_df =      '%s/faa_files_Kelp_and_HGT_B_GeneNumber.txt'                % wd
Kelp_ko_B_pct_df =      '%s/faa_files_Kelp_and_HGT_B_GeneNumber_pct.txt'            % wd
Tara_ko_B_num_df =      '%s/faa_files_Tara_NM_and_HGT_B_GeneNumber.txt'             % wd
Tara_ko_B_pct_df =      '%s/faa_files_Tara_NM_and_HGT_B_GeneNumber_pct.txt'         % wd
Kelp_ko_C_num_df =      '%s/faa_files_Kelp_and_HGT_C_GeneNumber.txt'                % wd
Kelp_ko_C_pct_df =      '%s/faa_files_Kelp_and_HGT_C_GeneNumber_pct.txt'            % wd
Tara_ko_C_num_df =      '%s/faa_files_Tara_NM_and_HGT_C_GeneNumber.txt'             % wd
Tara_ko_C_pct_df =      '%s/faa_files_Tara_NM_and_HGT_C_GeneNumber_pct.txt'         % wd

# output files
Kelp_ko_B_num_df_t =    '%s/faa_files_Kelp_and_HGT_B_GeneNumber_t.txt'              % wd
Kelp_ko_B_pct_df_t =    '%s/faa_files_Kelp_and_HGT_B_GeneNumber_pct_t.txt'          % wd
Tara_ko_B_num_df_t =    '%s/faa_files_Tara_NM_and_HGT_B_GeneNumber_t.txt'           % wd
Tara_ko_B_pct_df_t =    '%s/faa_files_Tara_NM_and_HGT_B_GeneNumber_pct_t.txt'       % wd
Kelp_ko_C_num_df_t =    '%s/faa_files_Kelp_and_HGT_C_GeneNumber_t.txt'              % wd
Kelp_ko_C_pct_df_t =    '%s/faa_files_Kelp_and_HGT_C_GeneNumber_pct_t.txt'          % wd
Tara_ko_C_num_df_t =    '%s/faa_files_Tara_NM_and_HGT_C_GeneNumber_t.txt'           % wd
Tara_ko_C_pct_df_t =    '%s/faa_files_Tara_NM_and_HGT_C_GeneNumber_pct_t.txt'       % wd

Kelp_ko_B_suppl_txt =   '%s/Kelp_ko_B_suppl.txt'                                    % wd
Tara_ko_B_suppl_txt =   '%s/Tara_ko_B_suppl.txt'                                    % wd
Kelp_ko_C_suppl_txt =   '%s/Kelp_ko_C_suppl.txt'                                    % wd
Tara_ko_C_suppl_txt =   '%s/Tara_ko_C_suppl.txt'                                    % wd


########################################################################################################################

# db files
ko00001_keg =                       '/Users/songweizhi/DB/KEGG_DB_important/ko00001.keg'

# read in db files
KO_description_A_dict, KO_description_B_dict, KO_description_C_dict, KO_description_D_dict, D2ABCD_dict, ko_list_C = get_ko2description_dict(ko00001_keg)

# transpose df
transpose_csv(Kelp_ko_B_num_df, Kelp_ko_B_num_df_t, '\t', 0, 0)
transpose_csv(Kelp_ko_B_pct_df, Kelp_ko_B_pct_df_t, '\t', 0, 0)
transpose_csv(Tara_ko_B_num_df, Tara_ko_B_num_df_t, '\t', 0, 0)
transpose_csv(Tara_ko_B_pct_df, Tara_ko_B_pct_df_t, '\t', 0, 0)
transpose_csv(Kelp_ko_C_num_df, Kelp_ko_C_num_df_t, '\t', 0, 0)
transpose_csv(Kelp_ko_C_pct_df, Kelp_ko_C_pct_df_t, '\t', 0, 0)
transpose_csv(Tara_ko_C_num_df, Tara_ko_C_num_df_t, '\t', 0, 0)
transpose_csv(Tara_ko_C_pct_df, Tara_ko_C_pct_df_t, '\t', 0, 0)

# get suppl table
get_suppl_table(Kelp_ko_B_num_df_t, Kelp_ko_B_pct_df_t, KO_description_B_dict, Kelp_ko_B_suppl_txt, 'B')
get_suppl_table(Tara_ko_B_num_df_t, Tara_ko_B_pct_df_t, KO_description_B_dict, Tara_ko_B_suppl_txt, 'B')
get_suppl_table(Kelp_ko_C_num_df_t, Kelp_ko_C_pct_df_t, KO_description_C_dict, Kelp_ko_C_suppl_txt, 'C')
get_suppl_table(Tara_ko_C_num_df_t, Tara_ko_C_pct_df_t, KO_description_C_dict, Tara_ko_C_suppl_txt, 'C')

# remove tmp files.
os.system('rm %s' % Kelp_ko_B_num_df_t)
os.system('rm %s' % Kelp_ko_B_pct_df_t)
os.system('rm %s' % Tara_ko_B_num_df_t)
os.system('rm %s' % Tara_ko_B_pct_df_t)
os.system('rm %s' % Kelp_ko_C_num_df_t)
os.system('rm %s' % Kelp_ko_C_pct_df_t)
os.system('rm %s' % Tara_ko_C_num_df_t)
os.system('rm %s' % Tara_ko_C_pct_df_t)

