import os
import glob
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
import s0_Kelp_bins_config
import seaborn as sns
from s0_suspicious_HGTs import suspicious_HGTs
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import s0_Kelp_bins_config as KelpCfg
from s0_suspicious_HGTs import suspicious_HGTs


############################################### input file and parameters ##############################################

genome_id_file = '/Users/songweizhi/Desktop/KelpBins/genome_id_list.txt'

day_list = ['050416', '060416', '080416',  # 2016-04: 05, 06, 08
             '070616', '140616',  # 2016-06: 07, 14
             '050716', '070716',  # 2016-07: 05, 07
             '110816', '170816',  # 2016-08: 11, 17
             '061016',  # 2016-10: 06
             '141216', '161216',  # 2016-12: 14, 16
             '070217', '080217',  # 2017-02: 07, 08
             '040417', '050417',  # 2017-04: 04, 05
             '050617', '130617']            # 2017-06: 05, 13

month_list = ['0416', '0616', '0716', '0816', '1016', '1216', '0217', '0417', '0617']
location_list = ['BH', 'BI', 'CB', 'SH']
host = 'ER'


wd = '/Users/songweizhi/Desktop/KelpBins/combined_pcofg'
grouping_file_folder =                  'GoodBins_0.5_0.05_pcofg_grouping'
PG_output_folder =                      'TT_90MGs_PG'
pwd_candidates_seq_file =               'GoodBins_0.5_0.05_all_combined_ffn.fasta'
output_prefix =                         ''
detection_ranks =                       'pcofg'
grouping_file_highest_rank =            '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_pcofg_grouping/GoodBins_0.5_0.05_p10_grouping.txt'
pwd_candidates_file_PG_normal_txt =     'GoodBins_0.5_0.05_PG_pcofg_normal.txt'
plot_all_genomes =                      False
pwd_group_to_taxon_file =               '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_pcofg_grouping/GoodBins_0.5_0.05_c15_group_to_taxon.txt'
grouping_file =                         '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_pcofg_grouping/GoodBins_0.5_0.05_c15_grouping.txt'

os.chdir(wd)

########################################################################################################################


# get the number of HGT at each month
sample_month_HGT_num = {}
for each in open(pwd_candidates_file_PG_normal_txt):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        concatenated_genes = '%s___%s' % (gene_1, gene_2)

        if concatenated_genes not in suspicious_HGTs:
            gene_1_genome = '_'.join(gene_1.split('_')[:-1])
            gene_2_genome = '_'.join(gene_2.split('_')[:-1])
            direction = each_split[6]

            donor_genome = direction.split('-->')[0]
            recipient_genome = direction.split('-->')[1]
            if '%)' in recipient_genome:
                recipient_genome = recipient_genome.split('(')[0]

            sample_date = recipient_genome.split('_')[2]
            sample_month = sample_date[2:]

            # get the number of HGT at each month
            if sample_month not in sample_month_HGT_num:
                sample_month_HGT_num[sample_month] = 1
            else:
                sample_month_HGT_num[sample_month] += 1


# get genome size dict
genome_size_dict = {}
for each_genome_size in open(KelpCfg.genome_size_file):
    if not each_genome_size.startswith('Genome'):
        each_genome_size_split = each_genome_size.strip().split('\t')
        genome_id = each_genome_size_split[0]
        genome_size = float(each_genome_size_split[1])
        genome_size_dict[genome_id] = genome_size


# get the total length of genome bins at each month
sample_month_genome_size_dict = {}
for genome in genome_size_dict:
    genome_size = genome_size_dict[genome]
    sample_month = genome.split('_')[2][2:]
    if sample_month not in sample_month_genome_size_dict:
        sample_month_genome_size_dict[sample_month] = genome_size
    else:
        sample_month_genome_size_dict[sample_month] += genome_size


print('The number of HGT per Mbp sequences at a monthly basis for all genome bins:')
for mth in month_list:
    mth_HGT_num = sample_month_HGT_num[mth]
    mth_genome_size = sample_month_genome_size_dict[mth]
    print('%s\t%s' % (mth, float("{0:.2f}".format(mth_HGT_num/mth_genome_size))))








