import os
import numpy as np
import pandas as pd
from scipy import stats
import Kelp_HGT_config


# test normality in python (scipy.stats.normaltest)
# test for equal variances: scipy.stats.levene and


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


##################################################### Prepare files ####################################################

# # transpose csv file
# transpose_csv(Kelp_HGT_config.Kelp_dbCAN_df_txt, Kelp_HGT_config.Kelp_dbCAN_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Tara_dbCAN_df_txt, Kelp_HGT_config.Tara_dbCAN_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Kelp_COG_df_txt, Kelp_HGT_config.Kelp_COG_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Tara_COG_df_txt, Kelp_HGT_config.Tara_COG_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Kelp_KEGG_df_txt, Kelp_HGT_config.Kelp_KEGG_df_txt_t, '\t', 0, 0)
# transpose_csv(Kelp_HGT_config.Tara_KEGG_df_txt, Kelp_HGT_config.Tara_KEGG_df_txt_t, '\t', 0, 0)


######################################################### dbCAN ########################################################

# get boxplot
# os.system('Rscript %s -i %s -o %s' % (Kelp_HGT_config.Boxplot_last1row, Kelp_HGT_config.Kelp_dbCAN_df_txt, Kelp_HGT_config.Kelp_dbCAN_df_png))
# os.system('Rscript %s -i %s -o %s' % (Kelp_HGT_config.Boxplot_last1row, Kelp_HGT_config.Tara_dbCAN_df_txt, Kelp_HGT_config.Tara_dbCAN_df_png))

