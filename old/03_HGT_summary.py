import statistics
import numpy as np
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import shutil


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def plot_donor_to_recipient_HGT_num(lol_in, sponge_type_list_sorted, plot_out):

    # turn list of list into arrary
    donor_to_recipient_HGT_num_arrary = np.array(lol_in)

    fig, ax = plt.subplots()
    im = ax.imshow(donor_to_recipient_HGT_num_arrary, interpolation='nearest')

    # show all ticks
    ax.set_xticks(np.arange(len(sponge_type_list_sorted)))
    ax.set_yticks(np.arange(len(sponge_type_list_sorted)))

    # label ticks
    ax.set_xticklabels(sponge_type_list_sorted, fontsize=7, fontstyle='italic')
    ax.set_yticklabels(sponge_type_list_sorted, fontsize=7, fontstyle='italic')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(sponge_type_list_sorted)):
        for j in range(len(sponge_type_list_sorted)):
            text = ax.text(j, i, donor_to_recipient_HGT_num_arrary[i, j],
                           ha="center", va="center", color="w", fontsize=8)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)

    # save the plot
    # ax.set_title("Harvest of local farmers (in tons/year)")
    fig.tight_layout()
    plt.savefig(plot_out, bbox_inches='tight', dpi=600)
    plt.close()
    plt.clf()


def box_plotter(num_list_in, label_list, plot_title, x_axis_label, y_axis_label, output_plot):

    # turn num list into arrary
    MAG_HGT_num_lol_arrary = [np.array(i) for i in num_list_in]

    # get plot
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(MAG_HGT_num_lol_arrary)

    # set x tick labels
    if len(label_list) <= 5:
        rotation_value = 0
    else:
        rotation_value = 270
    ax.set_xticklabels(label_list, rotation=rotation_value, fontsize=8)


    # set title, x and y label
    plt.title(plot_title)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='+', color='black', alpha=0.7, markersize=3)

    plt.tight_layout()
    fig.savefig(output_plot, bbox_inches='tight', dpi=300)
    plt.close()


def donor_to_recipient_HGT_heatmap(sponge_type_list_sorted, donor_to_recipient_HGT_num_matrix_file_absolute, donor_to_recipient_HGT_num_matrix_file_normalized, direction_by_host_count_dict, sponge_to_bin_total_size_dict, donor_to_recipient_HGT_num_plot_absolute, donor_to_recipient_HGT_num_plot_normalized, skip_plot):
    donor_to_recipient_HGT_num_matrix_file_absolute_handle = open(donor_to_recipient_HGT_num_matrix_file_absolute, 'w')
    donor_to_recipient_HGT_num_matrix_file_normalized_handle = open(donor_to_recipient_HGT_num_matrix_file_normalized, 'w')
    donor_to_recipient_HGT_num_matrix_file_absolute_handle.write(',%s\n' % ','.join(sponge_type_list_sorted))
    donor_to_recipient_HGT_num_matrix_file_normalized_handle.write(',%s\n' % ','.join(sponge_type_list_sorted))

    donor_to_recipient_HGT_num_absolute_lol = []
    donor_to_recipient_HGT_num_normalized_lol = []

    within_sponge_values = []
    between_sponge_values = []
    for sponge_donor in sponge_type_list_sorted:

        curret_donor_to_recipient_list_absolute = []
        curret_donor_to_recipient_list_normalized = []
        for sponge_recipient in sponge_type_list_sorted:
            donor_to_recipient = '%s-->%s' % (sponge_donor, sponge_recipient)

            # get donor_to_recipient_HGT_num
            if donor_to_recipient in direction_by_host_count_dict:
                donor_to_recipient_HGT_num_absolute = direction_by_host_count_dict[donor_to_recipient]
                donor_to_recipient_HGT_num_normalized = float("{0:.2f}".format(
                    donor_to_recipient_HGT_num_absolute / sponge_to_bin_total_size_dict[sponge_recipient]))
            else:
                donor_to_recipient_HGT_num_absolute = 0
                donor_to_recipient_HGT_num_normalized = 0

            if sponge_donor == sponge_recipient:
                within_sponge_values.append(donor_to_recipient_HGT_num_normalized)
            else:
                between_sponge_values.append(donor_to_recipient_HGT_num_normalized)

            curret_donor_to_recipient_list_absolute.append(donor_to_recipient_HGT_num_absolute)
            curret_donor_to_recipient_list_normalized.append(donor_to_recipient_HGT_num_normalized)

        donor_to_recipient_HGT_num_absolute_lol.append(curret_donor_to_recipient_list_absolute)
        donor_to_recipient_HGT_num_normalized_lol.append(curret_donor_to_recipient_list_normalized)

        current_donor_to_recipient_list_as_string_absolute   = [str(i) for i in curret_donor_to_recipient_list_absolute]
        current_donor_to_recipient_list_as_string_normalized = [str(i) for i in curret_donor_to_recipient_list_normalized]

        donor_to_recipient_HGT_num_matrix_file_absolute_handle.write('%s,%s\n' % (sponge_donor, ','.join(current_donor_to_recipient_list_as_string_absolute)))
        donor_to_recipient_HGT_num_matrix_file_normalized_handle.write('%s,%s\n' % (sponge_donor, ','.join(current_donor_to_recipient_list_as_string_normalized)))

    donor_to_recipient_HGT_num_matrix_file_absolute_handle.close()
    donor_to_recipient_HGT_num_matrix_file_normalized_handle.close()

    # get the plot
    if skip_plot is False:
        plot_donor_to_recipient_HGT_num(donor_to_recipient_HGT_num_absolute_lol, sponge_type_list_sorted, donor_to_recipient_HGT_num_plot_absolute)
        plot_donor_to_recipient_HGT_num(donor_to_recipient_HGT_num_normalized_lol, sponge_type_list_sorted, donor_to_recipient_HGT_num_plot_normalized)

    return within_sponge_values, between_sponge_values


def T_test(num_list_1, num_list_2, num_list_1_name, num_list_2_name):

    # turn list to arrary
    num_list_1_arrary = np.array(num_list_1)
    num_list_2_arrary = np.array(num_list_2)

    # get mean and stdev
    num_list_1_mean = statistics.mean(num_list_1_arrary)
    num_list_2_mean = statistics.mean(num_list_2_arrary)
    num_list_1_stdev = statistics.stdev(num_list_1_arrary)
    num_list_2_stdev = statistics.stdev(num_list_2_arrary)

    # perform t_test
    t_test= stats.ttest_ind(num_list_1_arrary,num_list_2_arrary)

    # turn num list to str list
    num_list_1_str = [str(i) for i in num_list_1]
    num_list_2_str = [str(i) for i in num_list_2]

    # report
    print('%s\tmean:%s\tstdev:%s' % (num_list_1_name, float("{0:.2f}".format(num_list_1_mean)), float("{0:.2f}".format(num_list_1_stdev))))
    print('%s\tmean:%s\tstdev:%s' % (num_list_2_name, float("{0:.2f}".format(num_list_2_mean)), float("{0:.2f}".format(num_list_2_stdev))))
    print('P-value: %s' % float("{0:.3f}".format(t_test.pvalue)))


####################################################### file I/O #######################################################

wd = '/Users/songweizhi/Desktop/SpongeEMP'

# file in
name_correlation_file =                                     '%s/0_file_in/bin_name_correlation.csv'                                         % wd
metadata_file =                                             '%s/0_file_in/SpongeEMP_metadata.txt'                                           % wd
genome_size_file =                                          '%s/0_file_in/SpongeEMP_all_genome_size.txt'                                    % wd
detected_HGTs_pcofg =                                       '%s/0_file_in/SpongeEMP_pcofg_detected_HGTs.txt'                                % wd
# detected_HGTs_pcofg =                                       '%s/0_file_in/SpongeEMP_pcofg_detected_HGTs_min_level_num_2.txt'                % wd
sponge_with_high_HGT_preferences =                          '%s/0_file_in/sponge_with_bin_size_top_7.txt'                                   % wd
skip_plot = True  # True or False

# file out
donor_to_recipient_HGT_num_matrix_file_absolute =           '%s/SpongeEMP_pcofg_donor_to_recipient_HGT_num_matrix_absolute.txt'             % wd
donor_to_recipient_HGT_num_matrix_file_absolute_high =      '%s/SpongeEMP_pcofg_donor_to_recipient_HGT_num_matrix_absolute_high.txt'        % wd
donor_to_recipient_HGT_num_matrix_file_normalized =         '%s/SpongeEMP_pcofg_donor_to_recipient_HGT_num_matrix_normalized.txt'           % wd
donor_to_recipient_HGT_num_matrix_file_normalized_high =    '%s/SpongeEMP_pcofg_donor_to_recipient_HGT_num_matrix_normalized_high.txt'      % wd
donor_to_recipient_HGT_num_plot_absolute =                  '%s/SpongeEMP_pcofg_donor_to_recipient_HGT_num_absolute.png'                    % wd
donor_to_recipient_HGT_num_plot_absolute_high =             '%s/SpongeEMP_pcofg_donor_to_recipient_HGT_num_absolute_high.png'               % wd
donor_to_recipient_HGT_num_plot_normalized =                '%s/SpongeEMP_pcofg_donor_to_recipient_HGT_num_normalized.png'                  % wd
donor_to_recipient_HGT_num_plot_normalized_high =           '%s/SpongeEMP_pcofg_donor_to_recipient_HGT_num_normalized_high.png'             % wd
plot_sponge_to_bin_num_and_size =                           '%s/sponge_to_bin_num_and_size.png'                                             % wd
MAG_to_HGT_num_file =                                       '%s/MAG_to_HGT_num.txt'                                                         % wd
MAG_to_HGT_num_plot =                                       '%s/MAG_to_HGT_num.png'                                                         % wd
host_to_HGT_num_file =                                      '%s/Host_to_HGT_num.txt'                                                        % wd


Donor_to_recipient_HGT_identity_plot_folder =               '%s/Donor_to_recipient_HGT_identity_plot'                                       % wd

if skip_plot is not True:
    force_create_folder(Donor_to_recipient_HGT_identity_plot_folder)


########################################################################################################################


# read in name correlation
new2raw_name_dict = {}
raw2new_name_dict = {}
for bin_name in open(name_correlation_file):
    raw_name = bin_name.strip().split(',')[0]
    new_name = bin_name.strip().split(',')[1]
    new2raw_name_dict[new_name] = raw_name
    raw2new_name_dict[raw_name] = new_name


# read in metadata
bin_to_source_dict = {}
bin_to_source2_dict = {}
sponge_to_bin_dict = {}
sponge_type_list = []
sponge_bin_list = []
seawater_bin_list = []
for each_bin in open(metadata_file):
    if not each_bin.startswith('Source'):
        each_bin_split = each_bin.strip().split(',')
        Source = each_bin_split[0]
        Source_2 = each_bin_split[1]
        MAG_ID = each_bin_split[5]
        bin_to_source_dict[MAG_ID] = Source

        # get sponge_type_list
        if (Source_2 not in sponge_type_list) and (Source_2 != 'TBA'):
            sponge_type_list.append(Source_2)

        if Source == 'seawater':
            bin_to_source2_dict[MAG_ID] = 'seawater'
            seawater_bin_list.append(MAG_ID)
        else:
            bin_to_source2_dict[MAG_ID] = Source_2
            sponge_bin_list.append(MAG_ID)

        # get sponge_to_bin_dict
        if Source == 'sponge':
            if Source_2 not in sponge_to_bin_dict:
                sponge_to_bin_dict[Source_2] = [MAG_ID]
            else:
                sponge_to_bin_dict[Source_2].append(MAG_ID)


# get bin size dict
bin_size_dict = {}
for genome in open(genome_size_file):
    if not genome.startswith('Genome	Size(Mbp)'):
        genome_sizes_plit = genome.strip().split('\t')
        genome_name = '.'.join(genome_sizes_plit[0].split('.')[:-1])
        genome_size = float(genome_sizes_plit[1])
        bin_size_dict[genome_name] = genome_size


sponge_to_bin_total_size_dict = {}
for sponge in sponge_to_bin_dict:
    current_sponge_bin_total_size = 0
    for current_sponge_bin in sponge_to_bin_dict[sponge]:
        current_sponge_bin_new_name = raw2new_name_dict[current_sponge_bin]
        current_sponge_bin_size = bin_size_dict[current_sponge_bin_new_name]
        current_sponge_bin_total_size += current_sponge_bin_size
    sponge_to_bin_total_size_dict[sponge] = current_sponge_bin_total_size


# get the total size of sponge and seawater bins
sponge_bins_size = 0
seawater_bins_size = 0
all_bins_size = 0
for genome_bin in bin_size_dict:
    genome_bin_size = bin_size_dict[genome_bin]
    genome_bin_raw_name = new2raw_name_dict[genome_bin]
    genome_bin_source = bin_to_source_dict[genome_bin_raw_name]
    all_bins_size += genome_bin_size
    if genome_bin_source == 'sponge':
        sponge_bins_size += genome_bin_size
    if genome_bin_source == 'seawater':
        seawater_bins_size += genome_bin_size


# get HGT summary
sponge_bins_HGT_num = 0
seawater_bins_HGT_num = 0
direction_by_host_count_dict = {}
bin_to_HGT_num_dict = {}
sponge_HGT_num_dict = {}
direction_by_host_HGT_identities_dict = {}
for each in open(detected_HGTs_pcofg):

    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        identity = float(each_split[2])
        direction = each_split[6]

        # get recipient_genome
        recipient_genome = direction.split('-->')[1]
        if '%)' in recipient_genome:
            recipient_genome = recipient_genome.split('(')[0]

        # get recipient_gene
        if recipient_genome == gene_1_genome:
            recipient_gene = gene_1
            donor_gene = gene_2
            donor_genome = gene_2_genome
        else:
            recipient_gene = gene_2
            donor_gene = gene_1
            donor_genome = gene_1_genome

        # get bin_to_HGT_num_dict
        if recipient_genome not in bin_to_HGT_num_dict:
            bin_to_HGT_num_dict[recipient_genome] = 1
        else:
            bin_to_HGT_num_dict[recipient_genome] += 1

        # get genome raw name
        recipient_genome_raw_name = new2raw_name_dict[recipient_genome]
        donor_genome_raw_name = new2raw_name_dict[donor_genome]

        # get recipient genome source
        recipient_genome_source     = bin_to_source_dict[recipient_genome_raw_name]
        recipient_genome_source2    = bin_to_source2_dict[recipient_genome_raw_name]
        donor_genome_source         = bin_to_source_dict[donor_genome_raw_name]
        donor_genome_source2        = bin_to_source2_dict[donor_genome_raw_name]

        # get the number of HGT for sponge and seawater bins
        if recipient_genome_source == 'sponge':
            sponge_bins_HGT_num += 1

            if recipient_genome_source2 not in sponge_HGT_num_dict:
                sponge_HGT_num_dict[recipient_genome_source2] = 1
            else:
                sponge_HGT_num_dict[recipient_genome_source2] += 1

        if recipient_genome_source == 'seawater':
            seawater_bins_HGT_num += 1

            if 'seawater' not in sponge_HGT_num_dict:
                sponge_HGT_num_dict['seawater'] = 1
            else:
                sponge_HGT_num_dict['seawater'] += 1

        direction_by_host = '%s-->%s' % (donor_genome_source2, recipient_genome_source2)

        # get direction_by_host_HGT_identities_dict
        if direction_by_host not in direction_by_host_HGT_identities_dict:
            direction_by_host_HGT_identities_dict[direction_by_host] = [identity]
        else:
            direction_by_host_HGT_identities_dict[direction_by_host].append(identity)


        if (recipient_genome_source == 'sponge') and (donor_genome_source == 'sponge'):

            if direction_by_host not in direction_by_host_count_dict:
                direction_by_host_count_dict[direction_by_host] = 1
            else:
                direction_by_host_count_dict[direction_by_host] += 1


######################################## get donor_to_recipient_HGT_num_matrix ########################################

# get sponge_with_high_HGT_preferences_list
sponge_with_high_HGT_preferences_list = []
for sponge in open(sponge_with_high_HGT_preferences):
    sponge_with_high_HGT_preferences_list.append(sponge.strip())

# sort sponge_type_list
sponge_type_list_sorted = sorted(sponge_type_list)
sponge_with_high_HGT_preferences_list_sorted = sorted(sponge_with_high_HGT_preferences_list)

within_sponge_values, between_sponge_values = donor_to_recipient_HGT_heatmap(sponge_type_list_sorted,
                                                                             donor_to_recipient_HGT_num_matrix_file_absolute,
                                                                             donor_to_recipient_HGT_num_matrix_file_normalized,
                                                                             direction_by_host_count_dict,
                                                                             sponge_to_bin_total_size_dict,
                                                                             donor_to_recipient_HGT_num_plot_absolute,
                                                                             donor_to_recipient_HGT_num_plot_normalized,
                                                                             skip_plot)

within_sponge_values_high, between_sponge_values_high = donor_to_recipient_HGT_heatmap(sponge_with_high_HGT_preferences_list_sorted,
                                                                                       donor_to_recipient_HGT_num_matrix_file_absolute_high,
                                                                                       donor_to_recipient_HGT_num_matrix_file_normalized_high,
                                                                                       direction_by_host_count_dict,sponge_to_bin_total_size_dict,
                                                                                       donor_to_recipient_HGT_num_plot_absolute_high,
                                                                                       donor_to_recipient_HGT_num_plot_normalized_high,
                                                                                       skip_plot)

################################ get the number and total size of bins for each sources ################################

source2_list_bin_num_no_seawater = [len(sponge_to_bin_dict[i]) for i in sponge_type_list_sorted]
source2_list_bin_total_size_no_seawater = [float("{0:.2f}".format(sponge_to_bin_total_size_dict[i])) for i in sponge_type_list_sorted]

source2_list_with_seawater = [i for i in sponge_type_list_sorted]
source2_list_with_seawater.append('seawater')

source2_list_bin_num_with_seawater = [i for i in source2_list_bin_num_no_seawater]
source2_list_bin_num_with_seawater.append(len(seawater_bin_list))


# get seawater_bins_total_size
seawater_bins_total_size = 0
for seawater_bins in seawater_bin_list:
    seawater_bins_size = bin_size_dict[seawater_bins]
    seawater_bins_total_size += seawater_bins_size
seawater_bins_total_size = float("{0:.2f}".format(seawater_bins_total_size))
source2_list_bin_total_size_with_seawater = [i for i in source2_list_bin_total_size_no_seawater]
source2_list_bin_total_size_with_seawater.append(seawater_bins_total_size)


# set width of bar
barWidth = 0.4

# set height of bar
bin_num_bar = source2_list_bin_num_with_seawater
bin_size_bar = source2_list_bin_total_size_with_seawater

# Set position of bar on X axis
r1 = np.arange(len(bin_num_bar))
r2 = [x + barWidth for x in r1]

if skip_plot is False:
    plt.bar(r1, bin_num_bar, color='orange', width=barWidth, edgecolor='white', label='Number')
    plt.bar(r2, bin_size_bar, color='lightgreen', width=barWidth, edgecolor='white', label='Size (Mbp)')
    plt.xticks([r + barWidth for r in range(len(bin_num_bar))], source2_list_with_seawater, rotation=270, fontsize=10)
    plt.title('The number and total size of bins')
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(plot_sponge_to_bin_num_and_size, dpi=300)
    plt.close()


####################################### export and plot MAG to HGT number stats ########################################

# export txt file
MAG_to_HGT_num_file_handle = open(MAG_to_HGT_num_file, 'w')
MAG_to_HGT_num_file_handle.write('Source,Bin,Num,Num_normalized(/Mbp)\n')
MAG_HGT_num_lol = []
for source2 in source2_list_with_seawater:

    if source2 != 'seawater':
        source2_bins = sponge_to_bin_dict[source2]
    else:
        source2_bins = seawater_bin_list

    current_MAG_HGT_num = []
    for MAG in source2_bins:
        MAG_new_name = raw2new_name_dict[MAG]
        if MAG_new_name in bin_to_HGT_num_dict:
            MAG_HGT_num = bin_to_HGT_num_dict[MAG_new_name]
            MAG_HGT_num_normalized = float("{0:.2f}".format(MAG_HGT_num/bin_size_dict[MAG_new_name]))
            current_MAG_HGT_num.append(MAG_HGT_num_normalized)
            MAG_to_HGT_num_file_handle.write('%s,%s,%s,%s\n' % (source2, MAG_new_name, MAG_HGT_num, MAG_HGT_num_normalized))

    MAG_HGT_num_lol.append(current_MAG_HGT_num)

MAG_to_HGT_num_file_handle.close()


# get box plot
if skip_plot is False:
    MAG_HGT_num_lol_arrary = [np.array(i) for i in MAG_HGT_num_lol]
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(MAG_HGT_num_lol_arrary)
    ax.set_xticklabels(source2_list_with_seawater, rotation=270, fontsize=8)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='+', color='black', alpha=0.7, markersize=3)

    plt.tight_layout()
    fig.savefig(MAG_to_HGT_num_plot, bbox_inches='tight', dpi=300)
    plt.close()


######################################## sponge bin size to HGT num correlation ########################################

# get sponge_HGT_num_list
host_to_HGT_num_file_handle = open(host_to_HGT_num_file, 'w')
host_to_HGT_num_file_handle.write('Host\tSize\tHGT_num/Mbp\n')
sponge_HGT_num_list = []
sponge_HGT_num_list_normalized = []
for source2 in source2_list_with_seawater:
    if source2 in sponge_HGT_num_dict:
        source2_HGT_num = sponge_HGT_num_dict[source2]

        if source2 != 'seawater':
            source2_HGT_num_normalized = source2_HGT_num/sponge_to_bin_total_size_dict[source2]
            host_to_HGT_num_file_handle.write('%s\t%s\t%s\n' % (source2, float("{0:.2f}".format(sponge_to_bin_total_size_dict[source2])), float("{0:.2f}".format(source2_HGT_num_normalized))))

        else:
            source2_HGT_num_normalized = source2_HGT_num/seawater_bins_total_size

    else:
        source2_HGT_num = 0
        source2_HGT_num_normalized = 0
        host_to_HGT_num_file_handle.write('%s\t%s\t%s\n' % (
        source2, float("{0:.2f}".format(sponge_to_bin_total_size_dict[source2])),float("{0:.2f}".format(source2_HGT_num_normalized))))

    sponge_HGT_num_list.append(source2_HGT_num)
    sponge_HGT_num_list_normalized.append(source2_HGT_num_normalized)
host_to_HGT_num_file_handle.close()

# # get the plot
# if skip_plot is False:
#     plt.scatter(np.array(source2_list_bin_total_size_with_seawater), np.array(sponge_HGT_num_list_normalized))
#     plt.show()
#     plt.savefig('/Users/songweizhi/Desktop/aaa.png', bbox_inches='tight', dpi=300)
#     plt.close()


########################################################################################################################

dornor_to_recipient_HGT_label_list = []
dornor_to_recipient_HGT_identity_lol = []
within_sponge_HGT_identity_list = []
between_sponge_HGT_identity_list = []
for sponge_recipient in sponge_with_high_HGT_preferences_list_sorted:

    # plot name
    current_recipient_dornor_to_recipient_HGT_identity_plot = '%s/HGT_identity_%s.png' % (Donor_to_recipient_HGT_identity_plot_folder, sponge_recipient)

    current_recipient_dornor_to_recipient_HGT_identity_lol = []
    for sponge_dornor in sponge_with_high_HGT_preferences_list_sorted:

        dornor_to_recipient = '%s-->%s' % (sponge_dornor, sponge_recipient)

        if dornor_to_recipient in direction_by_host_HGT_identities_dict:
            dornor_to_recipient_HGT_identity_list = direction_by_host_HGT_identities_dict[dornor_to_recipient]
        else:
            dornor_to_recipient_HGT_identity_list = []

        current_recipient_dornor_to_recipient_HGT_identity_lol.append(dornor_to_recipient_HGT_identity_list)

        dornor_to_recipient_HGT_label_list.append(dornor_to_recipient)
        dornor_to_recipient_HGT_identity_lol.append(dornor_to_recipient_HGT_identity_list)

        if sponge_dornor == sponge_recipient:
            within_sponge_HGT_identity_list.append(dornor_to_recipient_HGT_identity_list)
        else:
            between_sponge_HGT_identity_list.append(dornor_to_recipient_HGT_identity_list)

    # get plot
    if skip_plot is False:
        box_plotter(current_recipient_dornor_to_recipient_HGT_identity_lol, sponge_with_high_HGT_preferences_list_sorted, 'Recipient: %s' % sponge_recipient, 'Donor', 'Identity (%)', current_recipient_dornor_to_recipient_HGT_identity_plot)


################################# plot combined HGT identity for the top seven sponges #################################

within_sponge_HGT_identity_list_combined = []
for identity_list in within_sponge_HGT_identity_list:
    for identity in identity_list:
        within_sponge_HGT_identity_list_combined.append(identity)

between_sponge_HGT_identity_list_combined = []
for identity_list in between_sponge_HGT_identity_list:
    for identity in identity_list:
        between_sponge_HGT_identity_list_combined.append(identity)

if skip_plot is False:
    HGT_identity_all_sponges_combined = '%s/HGT_identity_all_sponges_combined.png' % Donor_to_recipient_HGT_identity_plot_folder
    box_plotter([within_sponge_HGT_identity_list_combined, between_sponge_HGT_identity_list_combined], ['Within sponge', 'Between sponges'], 'All sponges combined', '', 'Identity (%)', HGT_identity_all_sponges_combined)


#################################### get the number of HGTs detected from each host ####################################




######################################################## report ########################################################

print('\n========================== Report =============================\n')

print('The number of sponge types (without unknown sponge): %s' % (len(sponge_type_list) - 1))
print('The total size of all bins: %s Mbp'                      % float("{0:.2f}".format(all_bins_size)))
print()
print('Sponge bin number:             %s'                       % len(sponge_bin_list))
print('Sponge bin total size:         %s Mbp'                   % float("{0:.2f}".format(sponge_bins_size)))
print('Sponge bin #HGT:               %s'                       % sponge_bins_HGT_num)
print('Sponge bin #HGT/Mbp:           %s'                       % (float("{0:.2f}".format(sponge_bins_HGT_num/sponge_bins_size))))
print()
print('Seawater bin number:           %s'                       % len(seawater_bin_list))
print('Seawater bin total size:       %s Mbp'                   % float("{0:.2f}".format(seawater_bins_total_size)))
print('Seawater bin #HGT:             %s'                       % seawater_bins_HGT_num)
print('Seawater bin #HGT/Mbp:         %s'                       % (float("{0:.2f}".format(seawater_bins_HGT_num/seawater_bins_total_size))))

print('\nSummary of the number of HGT per Mbp sequences:')
T_test(within_sponge_values_high, between_sponge_values_high, 'Within sponge', 'Between sponges')

print('\nSummary of the identity of HGTs (Top seven sponges):')
T_test(within_sponge_HGT_identity_list_combined, between_sponge_HGT_identity_list_combined, 'Within sponge', 'Between sponges')

print('\n===============================================================\n')

