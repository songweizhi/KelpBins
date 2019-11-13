import statistics as stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import s0_Kelp_bins_config
from s0_suspicious_HGTs import suspicious_HGTs


def plot_list(HGT_per_Mbp_list, HGT_num_per_genome_boxplot):
    HGT_per_Mbp_list_with_class = []
    for each in HGT_per_Mbp_list:
        HGT_per_Mbp_list_with_class.append([each, 1])

    HGT_per_Mbp_list_with_class_df = pd.DataFrame(np.array(HGT_per_Mbp_list_with_class), columns=['num', 'label'])

    sns.set(style="whitegrid")  # whitegrid

    # get boxplot
    orient = 'v'
    ax = sns.boxplot(x="label", y="num", data=HGT_per_Mbp_list_with_class_df, orient=orient, width=0.5,
                     showfliers=False)
    ax = sns.stripplot(x="label", y="num", data=HGT_per_Mbp_list_with_class_df, orient=orient, color='orange', size=3,
                       jitter=0.23)
    ax.set(ylabel='Number of HGT/Mbp sequence', xlabel='Genome')

    # plt.show()
    plt.savefig(HGT_num_per_genome_boxplot, dpi=300)
    plt.close()


multi_level_detection = True
HGT_num_per_genome_boxplot = '/Users/songweizhi/Desktop/HGT_num_per_genome_boxplot.png'
HGT_num_per_genome_boxplot_all = '/Users/songweizhi/Desktop/HGT_num_per_genome_boxplot_all.png'


# store genome size in dict
genome_size_dict = {}
for each_size in open(s0_Kelp_bins_config.genome_size_file):
    if each_size.strip() != 'Genome\tSize(Mbp)':
        each_size_split = each_size.strip().split('\t')
        genome_file_name = each_size_split[0]
        genome_size_Mbp = float(each_size_split[1])
        genome_size_dict[genome_file_name] = genome_size_Mbp


# only focus on recipient
genome_recipient_num_dict = {}
for each_HGT in open(s0_Kelp_bins_config.HGT_PG_validated_txt_pcofg):

    if not each_HGT.startswith('Gene_1'):
        each_HGT_split = each_HGT.strip().split('\t')
        each_HGT_gene_1 = each_HGT_split[0]
        each_HGT_gene_2 = each_HGT_split[1]
        each_HGT_genome_1 = '_'.join(each_HGT_gene_1.split('_')[:-1])
        each_HGT_genome_2 = '_'.join(each_HGT_gene_2.split('_')[:-1])
        direction = each_HGT_split[5]
        concatenated_genes = '%s___%s' % (each_HGT_gene_1, each_HGT_gene_2)

        # only plot those with high confidence
        if concatenated_genes not in suspicious_HGTs:

            if multi_level_detection is True:
                direction = each_HGT_split[6]

            if '%)' in direction:
                direction = direction.split('(')[0]

            if each_HGT_genome_1 == direction.split('-->')[1]:
                recipient_genome = each_HGT_genome_1
                recipient_gene = each_HGT_gene_1
            else:
                recipient_gene = each_HGT_gene_2
                recipient_genome = each_HGT_genome_2


            if recipient_genome not in genome_recipient_num_dict:
                genome_recipient_num_dict[recipient_genome] = 1
            else:
                genome_recipient_num_dict[recipient_genome] += 1




for genome in genome_recipient_num_dict:
    genome_size = genome_size_dict[genome]
    HGT_num_per_Mbp = float("{0:.2f}".format(genome_recipient_num_dict[genome]/genome_size))

    print('%s\t%s\t%s' % (genome, genome_recipient_num_dict[genome], HGT_num_per_Mbp))




HGT_per_Mbp_list = []
HGT_per_Mbp_list_all_genome = []
for genome in genome_size_dict:
    genome_size = genome_size_dict[genome]
    HGT_num = 0
    if genome in genome_recipient_num_dict:
        HGT_num = genome_recipient_num_dict[genome]
    HGT_num_per_Mbp = float("{0:.2f}".format(HGT_num/genome_size))

    if HGT_num_per_Mbp > 0:
        HGT_per_Mbp_list.append(HGT_num_per_Mbp)

    HGT_per_Mbp_list_all_genome.append(HGT_num_per_Mbp)


print('The number of genomes with identified HGTs: %s' % len(genome_recipient_num_dict))

print(stats.median(HGT_per_Mbp_list))

plot_list(HGT_per_Mbp_list, HGT_num_per_genome_boxplot)
plot_list(HGT_per_Mbp_list_all_genome, HGT_num_per_genome_boxplot_all)



