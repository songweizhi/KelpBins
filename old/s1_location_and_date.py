import numpy as np
import matplotlib.pyplot as plt


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

##################################################### output file ######################################################

plot_out_day =   '/Users/songweizhi/Desktop/KelpBins/genome_per_sample_day.png'
plot_out_month = '/Users/songweizhi/Desktop/KelpBins/genome_per_sample_month.png'


########################################################################################################################

# get 528 genome is list
genome_528_list = []
for genome in open(genome_id_file):
    genome_528_list.append(genome.strip())


# get day matrix
genome_num_lol_day = []
for date in day_list:
    current_day_stats = []
    for location in location_list:
        tag = '%s_%s_%s' % (location, host, date)
        current_tag_genome_num = 0
        for genome in genome_528_list:
            if tag in genome:
                current_tag_genome_num += 1
        current_day_stats.append(current_tag_genome_num)
    current_date_stats_str = [str(i) for i in current_day_stats]
    genome_num_lol_day.append(current_day_stats)
    print('Day\t%s' % '\t'.join(current_date_stats_str))



genome_num_lol_day_array = np.array(genome_num_lol_day)


# plot
fig, ax = plt.subplots(figsize=(10,10))
im = ax.imshow(genome_num_lol_day_array, interpolation='none')

# show all ticks
ax.set_xticks(np.arange(len(location_list)))
ax.set_yticks(np.arange(len(day_list)))

# add labels
ax.set_xticklabels(location_list)
ax.set_yticklabels(day_list)

# Rotate x tick labels
plt.setp(ax.get_xticklabels(), rotation=0, ha="center", rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(day_list)):
    for j in range(len(location_list)):
        text = ax.text(j, i, genome_num_lol_day_array[i, j],
                       ha="center", va="center", color="w")

ax.set_title('Genome number per sample')
plt.savefig(plot_out_day, dpi=300)


########################################################################################################################

# get mth matrix
genome_num_lol_mth = []
for mth in month_list:
    current_day_stats = []
    for location in location_list:
        tag = '%s_%s_%s' % (location, host, mth)
        current_tag_genome_num = 0
        for genome in genome_528_list:
            genome_split = genome.split('_')
            genome_tag = '%s_%s_%s' % (genome_split[0], genome_split[1], genome_split[2][2:])
            if genome_tag == tag:
                current_tag_genome_num += 1
        current_day_stats.append(current_tag_genome_num)
    current_date_stats_str = [str(i) for i in current_day_stats]
    genome_num_lol_mth.append(current_day_stats)
    print('Mth\t%s' % '\t'.join(current_date_stats_str))


genome_num_lol_mth_array = np.array(genome_num_lol_mth)


# plot
fig, ax = plt.subplots(figsize=(7,7))
im = ax.imshow(genome_num_lol_mth_array, interpolation='none')

# show all ticks
ax.set_xticks(np.arange(len(location_list)))
ax.set_yticks(np.arange(len(month_list)))

# add labels
ax.set_xticklabels(location_list)
ax.set_yticklabels(month_list)

# Rotate x tick labels
plt.setp(ax.get_xticklabels(), rotation=0, ha="center", rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(month_list)):
    for j in range(len(location_list)):
        text = ax.text(j, i, genome_num_lol_mth_array[i, j],
                       ha="center", va="center", color="w")

ax.set_title('Genome number per sample')
plt.savefig(plot_out_month, dpi=300)
