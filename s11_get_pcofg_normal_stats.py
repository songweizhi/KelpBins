import numpy as np
import matplotlib.pyplot as plt

import s0_Kelp_bins_config


# store genome size in dict
genome_size_dict = {}
for each_size in open(s0_Kelp_bins_config.genome_size_file):
    if each_size.strip() != 'Genome\tSize(Mbp)':
        each_size_split = each_size.strip().split('\t')
        genome_file_name = each_size_split[0]
        genome_size_Mbp = float(each_size_split[1])
        genome_size_dict[genome_file_name] = genome_size_Mbp


genome_HGT_num_dict = {}
for each_HGT in open(s0_Kelp_bins_config.HGT_PG_validated_txt_pcofg):

    if not each_HGT.startswith('Gene_1'):
        each_HGT_split = each_HGT.strip().split('\t')
        each_HGT_gene_1 = each_HGT_split[0]
        each_HGT_gene_2 = each_HGT_split[1]
        each_HGT_genome_1 = '_'.join(each_HGT_gene_1.split('_')[:-1])
        each_HGT_genome_2 = '_'.join(each_HGT_gene_2.split('_')[:-1])

        if each_HGT_genome_1 not in genome_HGT_num_dict:
            genome_HGT_num_dict[each_HGT_genome_1] = 1
        else:
            genome_HGT_num_dict[each_HGT_genome_1] += 1

        if each_HGT_genome_2 not in genome_HGT_num_dict:
            genome_HGT_num_dict[each_HGT_genome_2] = 1
        else:
            genome_HGT_num_dict[each_HGT_genome_2] += 1


HGT_per_Mbp_list = []
for genome in genome_HGT_num_dict:
    genome_size = genome_size_dict[genome]
    HGT_num_per_Mbp = float("{0:.2f}".format(genome_HGT_num_dict[genome]/genome_size))
    for_print = '%s\t%s\t%s' % (genome, genome_HGT_num_dict[genome], HGT_num_per_Mbp)
    print(for_print)
    HGT_per_Mbp_list.append(HGT_num_per_Mbp)


print('The number of genomes with identified HGTs: %s' % len(genome_HGT_num_dict))


import seaborn as sns
import pandas as pd



print(HGT_per_Mbp_list)


HGT_per_Mbp_list_with_class = []

for each in HGT_per_Mbp_list:
    HGT_per_Mbp_list_with_class.append([each, 1])


print(HGT_per_Mbp_list_with_class)


HGT_per_Mbp_list_with_class_df = pd.DataFrame(np.array(HGT_per_Mbp_list_with_class), columns=['num', 'label'])
print(HGT_per_Mbp_list_with_class_df)



sns.set(style="whitegrid")  # whitegrid

# get boxplot
ax = sns.boxplot(x="num", y="label", data=HGT_per_Mbp_list_with_class_df, orient='h', width=0.5, showfliers=False)
ax = sns.stripplot(x="num", y="label", data=HGT_per_Mbp_list_with_class_df, color='orange', size=5, orient='h', jitter=0.23)
ax.set(xlabel='Number of HGT/Mbp sequence', ylabel='Genome')

plt.show()


# fig1, ax1 = plt.subplots()
# ax1.set_title('Basic Plot')
# ax1.boxplot(HGT_per_Mbp_list)
# plt.show()

