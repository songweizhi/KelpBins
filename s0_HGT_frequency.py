import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def plot_taxon_HGT_boxplot(identified_taxon_list_sorted, taxon_2_bin_HGT_num_normalized_dict, boxplot_png_filtered):

    # get lol for plot
    taxon_HGT_num_lol = []
    for identified_taxon in identified_taxon_list_sorted:
        taxon_HGT_num_lol.append(taxon_2_bin_HGT_num_normalized_dict[identified_taxon])

    # get boxplot
    MAG_HGT_num_lol_arrary = [np.array(i) for i in taxon_HGT_num_lol]
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(MAG_HGT_num_lol_arrary,
                    showfliers=False,
                    patch_artist=True,
                    whiskerprops=dict(color='lightblue', linewidth=2),
                    capprops=dict(color='lightblue'))

    # set the color pf box
    for box in bp['boxes']:
        box.set(linewidth=0)
        box.set_facecolor('lightblue')

    # for i in list(range(1, len(MAG_HGT_num_lol_arrary))):
    #     y = MAG_HGT_num_lol_arrary[i-1]
    #     x = np.random.normal(i, 0.08, len(y))
    #     plt.plot(x, y, '.', alpha=0.8, color='orange', markersize=5)

    ax.set_xticklabels(identified_taxon_list_sorted, fontsize=9, fontstyle='italic')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='+', color='black', alpha=0.7, markersize=3)

    ax.set(ylabel='#HGT/Mbp sequences')  # xlabel='Phylum'

    plt.tight_layout()
    fig.savefig(boxplot_png_filtered, bbox_inches='tight', dpi=300)
    plt.close()


######################################################## dRep99 ########################################################

wd = '/Users/songweizhi/Desktop/Kelp_NM'

# file in
dereplicated_bins =                     '%s/0_file_in/dereplicated_1030_bins.txt'                     % wd
genome_size_file =                      '%s/0_file_in/Kelp_NM_dRep99_all_genome_size.txt'             % wd
detected_HGTs_pcofg =                   '%s/0_file_in/Kelp_NM_dRep99_pcofg_detected_HGTs.txt'         % wd
GTDB_output_file_dRep =                 '%s/0_file_in/Kelp_NM_GTDB_r89.tsv'                           % wd

# parameters
circos_plot_min_HGT_num =               50
unclear_classification_level =          'f'
min_Bin =                               10
min_HGT =                               10
high_HGT_preference_cutoff =            20

keep_no_HGT_bins = False

# file out
pwd_matrix =                            '%s/SpongeEMP_circos_top_5.txt'             % wd
pwd_plot_circos =                       '%s/SpongeEMP_circos_top_5.png'             % wd
pwd_matrix_high_HGT_num =               '%s/SpongeEMP_circos_top_5_min_HGT_%s.txt'  % (wd, circos_plot_min_HGT_num)
pwd_plot_circos_high_HGT_num =          '%s/SpongeEMP_circos_top_5_min_HGT_%s.png'  % (wd, circos_plot_min_HGT_num)

rank_to_position_dict = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
taxon_rank_list = ['p']


########################################################################################################################

# read GTDB output into dict
taxon_assignment_dict = {}
for each_genome in open(GTDB_output_file_dRep):
    if not each_genome.startswith('user_genome'):
        each_split = each_genome.strip().split('\t')
        bin_name = each_split[0]
        assignment_full = []
        if len(each_split) == 1:
            assignment_full = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        elif (len(each_split) > 1) and (';' in each_split[1]):
            assignment = each_split[1].split(';')
            if len(assignment) == 7:
                assignment_full = assignment
            if len(assignment) == 6:
                assignment_full = assignment + ['s__']
            if len(assignment) == 5:
                assignment_full = assignment + ['g__', 's__']
            if len(assignment) == 4:
                assignment_full = assignment + ['f__', 'g__', 's__']
            if len(assignment) == 3:
                assignment_full = assignment + ['o__', 'f__', 'g__', 's__']
            if len(assignment) == 2:
                assignment_full = assignment + ['c__', 'o__', 'f__', 'g__', 's__']

        elif (len(each_split) > 1) and (';' not in each_split[1]):
            assignment_full = [each_split[1]] + ['p__', 'c__', 'o__', 'f__', 'g__', 's__']

        # store in dict
        taxon_assignment_dict[bin_name] = assignment_full


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
MAGs_with_high_HGT_preference = []
for genome_bin in bin_to_HGT_num_dict:
    bin_HGT_num = bin_to_HGT_num_dict[genome_bin]
    bin_size = bin_size_dict[genome_bin]
    bin_HGT_num_normalized = float("{0:.2f}".format(bin_HGT_num/bin_size))
    bin_to_HGT_num_normalized_dict[genome_bin] = bin_HGT_num_normalized

    if bin_HGT_num_normalized >= high_HGT_preference_cutoff:
        MAGs_with_high_HGT_preference.append(genome_bin)


for taxon_rank in taxon_rank_list:
    specified_rank_pos = rank_to_position_dict[taxon_rank]

    Kelp_taxon_HGT_freq_plot = '%s/Kelp_taxon_HGT_freq_%s_minBin_%s_minHGT_%s.png' % (wd, taxon_rank, min_Bin, min_HGT)
    Tara_taxon_HGT_freq_plot = '%s/Tara_taxon_HGT_freq_%s_minBin_%s_minHGT_%s.png' % (wd, taxon_rank, min_Bin, min_HGT)

    identified_Kelp_taxon_list = []
    identified_Tara_taxon_list = []
    taxon_2_Kelp_genome_dict = {}
    taxon_2_Tara_genome_dict = {}
    taxon_2_Kelp_HGT_num_dict = {}
    taxon_2_Tara_HGT_num_dict = {}
    taxon_2_Kelp_genome_HGT_freq_dict = {}
    taxon_2_Tara_genome_HGT_freq_dict = {}
    for each_TaxonAssign in taxon_assignment_dict:

        specified_rank_id = taxon_assignment_dict[each_TaxonAssign][specified_rank_pos]

        if '_Refined_' in each_TaxonAssign:  # for Kelp bins

            if specified_rank_id not in taxon_2_Kelp_genome_dict:
                taxon_2_Kelp_genome_dict[specified_rank_id] = [each_TaxonAssign]
                if each_TaxonAssign in bin_to_HGT_num_normalized_dict:
                    taxon_2_Kelp_genome_HGT_freq_dict[specified_rank_id] = [bin_to_HGT_num_normalized_dict[each_TaxonAssign]]
                    taxon_2_Kelp_HGT_num_dict[specified_rank_id] = bin_to_HGT_num_dict[each_TaxonAssign]
                else:
                    if keep_no_HGT_bins is True:
                        taxon_2_Kelp_genome_HGT_freq_dict[specified_rank_id] = [0]
                    else:
                        taxon_2_Kelp_genome_HGT_freq_dict[specified_rank_id] = []
                    taxon_2_Kelp_HGT_num_dict[specified_rank_id] = 0
            else:
                taxon_2_Kelp_genome_dict[specified_rank_id].append(each_TaxonAssign)
                if each_TaxonAssign in bin_to_HGT_num_normalized_dict:
                    taxon_2_Kelp_genome_HGT_freq_dict[specified_rank_id].append(bin_to_HGT_num_normalized_dict[each_TaxonAssign])
                    taxon_2_Kelp_HGT_num_dict[specified_rank_id] += bin_to_HGT_num_dict[each_TaxonAssign]

                else:
                    if keep_no_HGT_bins is True:
                        taxon_2_Kelp_genome_HGT_freq_dict[specified_rank_id].append(0)

            if specified_rank_id not in identified_Kelp_taxon_list:
                identified_Kelp_taxon_list.append(specified_rank_id)

        else:  # for Tara bins

            if specified_rank_id not in taxon_2_Tara_genome_dict:
                taxon_2_Tara_genome_dict[specified_rank_id] = [each_TaxonAssign]
                if each_TaxonAssign in bin_to_HGT_num_normalized_dict:
                    taxon_2_Tara_genome_HGT_freq_dict[specified_rank_id] = [bin_to_HGT_num_normalized_dict[each_TaxonAssign]]
                    taxon_2_Tara_HGT_num_dict[specified_rank_id] = bin_to_HGT_num_dict[each_TaxonAssign]
                else:
                    if keep_no_HGT_bins is True:
                        taxon_2_Tara_genome_HGT_freq_dict[specified_rank_id] = [0]
                    else:
                        taxon_2_Tara_genome_HGT_freq_dict[specified_rank_id] = []
                    taxon_2_Tara_HGT_num_dict[specified_rank_id] = 0
            else:
                taxon_2_Tara_genome_dict[specified_rank_id].append(each_TaxonAssign)
                if each_TaxonAssign in bin_to_HGT_num_normalized_dict:
                    taxon_2_Tara_genome_HGT_freq_dict[specified_rank_id].append(bin_to_HGT_num_normalized_dict[each_TaxonAssign])
                    taxon_2_Tara_HGT_num_dict[specified_rank_id] += bin_to_HGT_num_dict[each_TaxonAssign]
                else:
                    if keep_no_HGT_bins is True:
                        taxon_2_Tara_genome_HGT_freq_dict[specified_rank_id].append(0)

            if specified_rank_id not in identified_Tara_taxon_list:
                identified_Tara_taxon_list.append(specified_rank_id)

    # sort identified_taxon_list
    identified_Kelp_taxon_list_sorted = sorted(identified_Kelp_taxon_list)
    identified_Tara_taxon_list_sorted = sorted(identified_Tara_taxon_list)


    taxon_2_Kelp_genome_HGT_freq_dict_non_zero = {}
    for i in taxon_2_Kelp_genome_HGT_freq_dict:
        if (taxon_2_Kelp_genome_HGT_freq_dict[i] != []) and (len(taxon_2_Kelp_genome_dict[i]) >= min_Bin) and (taxon_2_Kelp_HGT_num_dict[i] >= min_HGT):
            taxon_2_Kelp_genome_HGT_freq_dict_non_zero[i] = taxon_2_Kelp_genome_HGT_freq_dict[i]

    taxon_2_Tara_genome_HGT_freq_dict_non_zero = {}
    for j in taxon_2_Tara_genome_HGT_freq_dict:
        if (taxon_2_Tara_genome_HGT_freq_dict[j] != []) and (len(taxon_2_Tara_genome_dict[j]) >= min_Bin) and (taxon_2_Tara_HGT_num_dict[j] >= min_HGT):
            taxon_2_Tara_genome_HGT_freq_dict_non_zero[j] = taxon_2_Tara_genome_HGT_freq_dict[j]

    qualified_Kelp_taxon_list_sorted = sorted([i for i in taxon_2_Kelp_genome_HGT_freq_dict_non_zero])
    qualified_Tara_taxon_list_sorted = sorted([i for i in taxon_2_Tara_genome_HGT_freq_dict_non_zero])

    plot_taxon_HGT_boxplot(qualified_Kelp_taxon_list_sorted, taxon_2_Kelp_genome_HGT_freq_dict_non_zero, Kelp_taxon_HGT_freq_plot)
    plot_taxon_HGT_boxplot(qualified_Tara_taxon_list_sorted, taxon_2_Tara_genome_HGT_freq_dict_non_zero, Tara_taxon_HGT_freq_plot)
