import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


######################################################## dRep99 ########################################################

wd = '/Users/songweizhi/Desktop/Kelp_NM'

# file in
dereplicated_bins =                     '%s/0_file_in/dereplicated_1030_bins.txt'                     % wd
genome_size_file =                      '%s/0_file_in/Kelp_NM_dRep99_all_genome_size.txt'             % wd
detected_HGTs_pcofg =                   '%s/0_file_in/Kelp_NM_dRep99_pcofg_detected_HGTs.txt'         % wd
mag_grouping =                          '%s/0_file_in/Kelp_NM_dRep99_p32_grouping.txt'                % wd

# file out
plot_out = '%s/HGT_freq.svg'                % wd

########################################################################################################################

# read GTDB output into dict
taxon_assignment_dict = {}
for each_genome in open(mag_grouping):
    each_genome_split = each_genome.strip().split(',')
    taxon_assignment_dict[each_genome_split[1]] = each_genome_split[2].split('__')[1]


# get bin size dict
bin_size_dict = {}
for genome in open(genome_size_file):
    if not genome.startswith('Genome	Size(Mbp)'):
        genome_sizes_plit = genome.strip().split('\t')
        genome_name = '.'.join(genome_sizes_plit[0].split('.')[:-1])
        genome_size = float(genome_sizes_plit[1])
        bin_size_dict[genome_name] = genome_size


# get genome_to_HGT_num_dict
bin_to_HGT_num_dict = {}
for each in open(detected_HGTs_pcofg):
    if not each.startswith('Gene_1'):
        direction = each.strip().split('\t')[6]

        # get recipient_genome
        recipient_genome = direction.split('-->')[1]
        if '%)' in recipient_genome:
            recipient_genome = recipient_genome.split('(')[0]

        if recipient_genome not in bin_to_HGT_num_dict:
            bin_to_HGT_num_dict[recipient_genome] = 1
        else:
            bin_to_HGT_num_dict[recipient_genome] += 1

bin_to_HGT_num_normalized_dict = {}
for genome_bin in bin_to_HGT_num_dict:
    bin_HGT_num = bin_to_HGT_num_dict[genome_bin]
    bin_size = bin_size_dict[genome_bin]
    bin_HGT_num_normalized = float("{0:.2f}".format(bin_HGT_num/bin_size))
    bin_to_HGT_num_normalized_dict[genome_bin] = bin_HGT_num_normalized

print(bin_to_HGT_num_normalized_dict)
print(len(bin_to_HGT_num_normalized_dict))


kelp_mag_HGT_num_list = []
tara_mag_HGT_num_list = []
for hgt_num in bin_to_HGT_num_normalized_dict:
    if '_Refined_' in hgt_num:
        kelp_mag_HGT_num_list.append(bin_to_HGT_num_normalized_dict[hgt_num])
    else:
        tara_mag_HGT_num_list.append(bin_to_HGT_num_normalized_dict[hgt_num])


print(len(kelp_mag_HGT_num_list))
print(len(tara_mag_HGT_num_list))
print(sum(kelp_mag_HGT_num_list)/len(kelp_mag_HGT_num_list))
print(sum(tara_mag_HGT_num_list)/len(tara_mag_HGT_num_list))

