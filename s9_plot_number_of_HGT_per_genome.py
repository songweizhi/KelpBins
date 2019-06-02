import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import s0_Kelp_bins_config as KelpCfg
from s0_suspicious_HGTs import suspicious_HGTs


genome_size_dict = {}
for each_genome_size in open(KelpCfg.genome_size_file):
    if not each_genome_size.startswith('Genome'):
        each_genome_size_split = each_genome_size.strip().split('\t')
        genome_id = each_genome_size_split[0]
        genome_size = float(each_genome_size_split[1])
        genome_size_dict[genome_id] = genome_size


genome_to_HGT_num_dict = {}
genome_set = set()
for hgt_pcofg in open(KelpCfg.HGT_PG_validated_txt_pcofg):
    if not hgt_pcofg.startswith('Gene_1'):
        hgt_pcofg_split = hgt_pcofg.strip().split('\t')
        gene_1 = hgt_pcofg_split[0]
        gene_2 = hgt_pcofg_split[1]
        concatenated_genes = '%s___%s' % (gene_1, gene_1)

        # only plot those with high confidence
        if concatenated_genes not in suspicious_HGTs:
            genome_1 = '_'.join(gene_1.split('_')[:-1])
            genome_2 = '_'.join(gene_2.split('_')[:-1])

            if genome_1 not in genome_to_HGT_num_dict:
                genome_to_HGT_num_dict[genome_1] = 1
            else:
                genome_to_HGT_num_dict[genome_1] += 1

            if genome_2 not in genome_to_HGT_num_dict:
                genome_to_HGT_num_dict[genome_2] = 1
            else:
                genome_to_HGT_num_dict[genome_2] += 1

            genome_set.add(genome_1)
            genome_set.add(genome_2)


genome_list = sorted([i for i in genome_set])
genome_HGT_num = [genome_to_HGT_num_dict[i] for i in genome_list]
genome_HGT_num_norm = [float("{0:.2f}".format((genome_to_HGT_num_dict[i]/genome_size_dict[i]))) for i in genome_list]


print(genome_list)
print(genome_HGT_num)
print(genome_HGT_num_norm)
xticks_fontsize_genome = 4
pwd_candidates_file_ET_validated_STAT_png = '/Users/songweizhi/Desktop/number_of_HGT_per_genome.png'
# set figure size
plt.figure(figsize=(20, 10))
x_range_genome = range(len(genome_HGT_num))

# subplot 1
plt.subplot(211)
plt.bar(x_range_genome, genome_HGT_num, alpha=0.5, linewidth=0, color='g', align='center')
plt.xticks(x_range_genome, genome_list, rotation=270, fontsize=xticks_fontsize_genome,
           horizontalalignment='center')
plt.axis('tight')
plt.ylabel('Number of HGT')


# subplot 3
plt.subplot(212)
plt.bar(x_range_genome, genome_HGT_num_norm, alpha=0.5, linewidth=0, color='g', align='center')
plt.xticks(x_range_genome, genome_list, rotation=270, fontsize=xticks_fontsize_genome,
           horizontalalignment='center')
plt.axis('tight')
plt.xlabel('Genome')
plt.ylabel('Number of HGT / Mbp sequences')


# plot layout
plt.subplots_adjust(wspace=0.2, top=0.95)

# save plot
plt.tight_layout()
plt.savefig(pwd_candidates_file_ET_validated_STAT_png, dpi=300)



