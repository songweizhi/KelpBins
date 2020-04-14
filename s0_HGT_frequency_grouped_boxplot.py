import pandas as pd
import seaborn as sns
import numpy as np
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
plot_out =                              '%s/HGT_freq.svg'                                             % wd

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


kelp_phylum_to_plot = ['Actinobacteriota', 'Bacteroidota', 'Planctomycetota', 'Proteobacteria']
Tara_phylum_to_plot = ['Actinobacteriota', 'Bacteroidota', 'Chloroflexota', 'Latescibacterota', 'Marinisomatota', 'Myxococcota', 'Planctomycetota', 'Proteobacteria', 'Thermoplasmatota', 'Verrucomicrobiota']

MAG_attributes_lol = []
for mag in bin_to_HGT_num_normalized_dict:
    mag_taxon = taxon_assignment_dict[mag]
    mag_source = 'Planktonic'
    if '_Refined_' in mag:
        mag_source = 'Kelp-associated'

    MAG_attribute_list = [mag, bin_to_HGT_num_normalized_dict[mag], mag_source, mag_taxon]

    if (mag_source == 'Kelp-associated') and (mag_taxon in kelp_phylum_to_plot):
        MAG_attributes_lol.append(MAG_attribute_list)
    if (mag_source == 'Planktonic') and (mag_taxon in Tara_phylum_to_plot):
        MAG_attributes_lol.append(MAG_attribute_list)


box_order =   ['Kelp-associated', 'Planktonic']
color_order = ['orange',          'lightblue']
MAG_attributes_df = pd.DataFrame(MAG_attributes_lol, columns=['MAG', 'Freq', 'Lifestyle', 'Phylum'])

chart = sns.boxplot(data=MAG_attributes_df, x="Phylum", y="Freq",order=Tara_phylum_to_plot,
            hue="Lifestyle", hue_order=box_order, palette=color_order,
            showfliers=False)

chart.set_xticklabels(chart.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.tight_layout()
plt.savefig(plot_out, bbox_inches='tight', dpi=300)
plt.close()
plt.clf()
