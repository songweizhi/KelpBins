import os
import numpy as np
import pandas as pd
from scipy import stats
import Kelp_HGT_config


Get_KEGG_boxplot = '''

# get plot
cd /Users/songweizhi/Desktop/Kelp_NM/KEGG
Rscript ~/PycharmProjects/BioSAK/BioSAK/KEGG_boxplot_last1row.R -i faa_files_Kelp_and_HGT_B_GeneNumber_pct.txt -o Kelp_and_HGT_KEGG_C_pct.png
Rscript ~/PycharmProjects/BioSAK/BioSAK/KEGG_boxplot_last1row.R -i faa_files_Tara_NM_and_HGT_B_GeneNumber_pct.txt -o Tara_and_HGT_KEGG_C_func_pct.png

'''
