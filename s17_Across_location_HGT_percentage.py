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



####################################################### file I/O #######################################################

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


########################################################################################################################

os.chdir(wd)


taxon_to_group_id_dict = {}
for group in open(pwd_group_to_taxon_file):
    group_id = group.strip().split(',')[0]
    group_taxon = group.strip().split(',')[1]
    taxon_to_group_id_dict[group_id] = group_taxon


# get genome to taxon dict
genome_to_taxon_dict = {}
for genome in open(grouping_file):
    group_id2 = genome.strip().split(',')[0]
    genome_name = genome.strip().split(',')[1]
    genome_to_taxon_dict[genome_name] = taxon_to_group_id_dict[group_id2]

within_group_HGT_num = 0
between_group_HGT_num = 0
taxon_in_between_transfer_dict = {}
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

            if (recipient_genome in genome_to_taxon_dict) and (donor_genome in genome_to_taxon_dict):
                recipient_taxon = genome_to_taxon_dict[recipient_genome]
                donor_taxon = genome_to_taxon_dict[donor_genome]

                if recipient_taxon not in taxon_in_between_transfer_dict:
                    if recipient_taxon == donor_taxon:
                        taxon_in_between_transfer_dict[recipient_taxon] = [1, 0]
                        within_group_HGT_num += 1
                    else:
                        taxon_in_between_transfer_dict[recipient_taxon] = [0, 1]
                        between_group_HGT_num += 1

                else:
                    if recipient_taxon == donor_taxon:
                        taxon_in_between_transfer_dict[recipient_taxon][0] += 1
                        within_group_HGT_num += 1
                    else:
                        taxon_in_between_transfer_dict[recipient_taxon][1] += 1
                        between_group_HGT_num += 1


# for each in taxon_in_between_transfer_dict:
#     print('%s\t%s\t%s' % (each, taxon_in_between_transfer_dict[each][0], taxon_in_between_transfer_dict[each][1]))

print('Within group HGT percentage: %s' % float("{0:.2f}".format(within_group_HGT_num*100/(within_group_HGT_num + between_group_HGT_num))) + '%')



####################################################### location #######################################################

within_location_HGT_num = 0
between_location_HGT_num = 0
taxon_in_between_transfer_dict = {}
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

            donor_location = donor_genome.split('_')[0]
            recipient_location = recipient_genome.split('_')[0]

            if recipient_location == donor_location:
                within_location_HGT_num += 1
            else:
                between_location_HGT_num += 1


print('Across location HGT percentage: %s' % float("{0:.2f}".format(within_location_HGT_num*100/(within_location_HGT_num + between_location_HGT_num))) + '%')


