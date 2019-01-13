import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


HGT_PG_validated_o = '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_o34_HGTs_PG_validated.txt'
pwd_iden_distribution_plot_PG = '/Users/songweizhi/Desktop/pwd_iden_distribution_plot_PG.png'


# get identity list
identity_list = []
for hgt_o in open(HGT_PG_validated_o):
    if not hgt_o.startswith('Gene_1'):
        hgt_o_split = hgt_o.strip().split('\t')
        identity = float(hgt_o_split[4])
        identity_list.append(identity)


# transfer identity to genetic variation
genetic_variation_list = []
for each in identity_list:
    genetic_variation_list.append(100 -each)


# Get plot
num_bins = 35
plt.hist(genetic_variation_list, num_bins, alpha=0.6, normed=0, linewidth=0, color='g', rwidth=0.95)
plt.title('Genetic variation of identified HGTs')
plt.xlabel('Genetic variation (%)')
plt.ylabel('Number of HGT')
plt.tight_layout()
plt.savefig(pwd_iden_distribution_plot_PG, bbox_inches='tight', dpi=300)
plt.close()
plt.clf()

