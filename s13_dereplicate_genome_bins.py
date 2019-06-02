import os
import shutil
import numpy as np
import matplotlib.pyplot as plt


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


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def get_GTDB_LSA(query_genome_list, combined_GTDB_summary_files, full_lineage):

    # get genome_to_taxon_dict
    genome_to_taxon_dict = {}
    for genome_taxon in open(combined_GTDB_summary_files):
        if not genome_taxon.startswith('user_genome'):
            genome_taxon_split = genome_taxon.strip().split('\t')
            genome_to_taxon_dict[genome_taxon_split[0]] = genome_taxon_split[1]

    # get last common ancestor
    taxon_rank_list = ['d', 'p', 'c', 'o', 'f', 'g', 's']

    LSA_determined = False
    last_common_ancestor = None
    for taxon_rank in taxon_rank_list[::-1]:
        current_rank_index = taxon_rank_list.index(taxon_rank)

        current_rank_member_list = []
        for query_genome in query_genome_list:
            query_genome_taxon = genome_to_taxon_dict[query_genome].split(';')
            current_rank_member_list.append(query_genome_taxon[current_rank_index])

        current_rank_member_list_uniq = unique_list_elements(current_rank_member_list)

        if (LSA_determined == False) and (len(current_rank_member_list_uniq) == 1) and (current_rank_member_list_uniq[0] != '%s__' % taxon_rank):
            last_common_ancestor = current_rank_member_list_uniq[0]
            LSA_determined = True

        elif (full_lineage == True) and (LSA_determined == True):
            last_common_ancestor = current_rank_member_list_uniq[0] + ';' + last_common_ancestor

        if (taxon_rank == 'd') and (len(current_rank_member_list_uniq) > 1):
            last_common_ancestor = 'life'

    return last_common_ancestor


###################################################### input files #####################################################

GTDB_results = '/Users/songweizhi/Desktop/KelpBins/gtdbtk.bac120.classification_528_r86.tsv'
drep_Cdb_file = '/Users/songweizhi/Desktop/KelpBins/Cdb.csv'
get_plot = True
genome_culster_summary_folder = '/Users/songweizhi/Desktop/genome_culster_summary_folder'


######################################################## config ########################################################

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


# define output file name
pwd_cluster_to_genome_member_file =     '%s/cluster_to_genome_member.txt'   % genome_culster_summary_folder
pwd_cluster_to_taxon_file =             '%s/cluster_to_taxon.txt'           % genome_culster_summary_folder
pwd_genome_num_per_sample_day_folder =  '%s/genome_num_per_sample_day'      % genome_culster_summary_folder
pwd_genome_num_per_sample_mth_folder =  '%s/genome_num_per_sample_mth'      % genome_culster_summary_folder


# create folder
force_create_folder(genome_culster_summary_folder)
os.mkdir(pwd_genome_num_per_sample_day_folder)
os.mkdir(pwd_genome_num_per_sample_mth_folder)


genome_to_cluster_dict = {}
cluster_to_genome_dict = {}
for genome_cluster in open(drep_Cdb_file):

    if not genome_cluster.startswith('primary_cluster'):
        genome_cluster_split = genome_cluster.strip().split(',')
        genome_id = genome_cluster_split[1]
        cluster_id = genome_cluster_split[0]

        # get genome_to_cluster_dict
        genome_to_cluster_dict[genome_id] = cluster_id

        # get cluster_to_genome_dict
        if cluster_id not in cluster_to_genome_dict:
            cluster_to_genome_dict[cluster_id] = [genome_id]
        else:
            cluster_to_genome_dict[cluster_id].append(genome_id)


cluster_to_genome_member_file_handle = open(pwd_cluster_to_genome_member_file, 'w')
cluster_to_taxon_file_handle = open(pwd_cluster_to_taxon_file, 'w')
for cluster in cluster_to_genome_dict:
    current_cluster_genome_list = cluster_to_genome_dict[cluster]
    current_cluster_genome_list_no_ext = ['.'.join(i.split('.')[:-1]) for i in current_cluster_genome_list]

    #################################### get LSA for each cluster according GTDB ###################################

    current_cluster_taxon = get_GTDB_LSA(current_cluster_genome_list_no_ext, GTDB_results, full_lineage=True)

    cluster_to_genome_member_file_handle.write('Cluster%s:\t%s\n' % (cluster, '\t'.join(current_cluster_genome_list)))
    cluster_to_taxon_file_handle.write('Cluster%s:\t%s\n' % (cluster, current_cluster_taxon))


    if (len(current_cluster_genome_list) > 1) and (get_plot == True):

        plot_out_day_current_cluster = '%s/Cluster%s_genome_num_per_sample_day.png' % (pwd_genome_num_per_sample_day_folder, cluster)
        plot_out_mth_current_cluster = '%s/Cluster%s_genome_num_per_sample_month.png' % (pwd_genome_num_per_sample_mth_folder, cluster)

        current_cluster_taxon_lowest = get_GTDB_LSA(current_cluster_genome_list_no_ext, GTDB_results, full_lineage=False)

        ###################################### plot genome number per sample (day) #####################################

        # get day matrix
        genome_num_lol_day = []
        for date in day_list:
            current_day_stats = []
            for location in location_list:
                tag = '%s_%s_%s' % (location, host, date)
                current_tag_genome_num = 0
                for genome in current_cluster_genome_list:
                    if tag in genome:
                        current_tag_genome_num += 1
                current_day_stats.append(current_tag_genome_num)
            current_date_stats_str = [str(i) for i in current_day_stats]
            genome_num_lol_day.append(current_day_stats)
            #print('Day\t%s' % '\t'.join(current_date_stats_str))

        genome_num_lol_day_array = np.array(genome_num_lol_day)

        # plot
        fig, ax = plt.subplots(figsize=(10, 10))
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

        ax.set_title(current_cluster_taxon_lowest)
        plt.savefig(plot_out_day_current_cluster, dpi=300)
        plt.close('all')

    ###################################### plot genome number per sample (mth) #####################################

        # get mth matrix
        genome_num_lol_mth = []
        for mth in month_list:
            current_day_stats = []
            for location in location_list:
                tag = '%s_%s_%s' % (location, host, mth)
                current_tag_genome_num = 0
                for genome in current_cluster_genome_list:
                    genome_split = genome.split('_')
                    genome_tag = '%s_%s_%s' % (genome_split[0], genome_split[1], genome_split[2][2:])
                    if genome_tag == tag:
                        current_tag_genome_num += 1
                current_day_stats.append(current_tag_genome_num)
            current_date_stats_str = [str(i) for i in current_day_stats]
            genome_num_lol_mth.append(current_day_stats)
            #print('Mth\t%s' % '\t'.join(current_date_stats_str))

        genome_num_lol_mth_array = np.array(genome_num_lol_mth)

        # plot
        fig, ax = plt.subplots(figsize=(7, 7))
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

        ax.set_title(current_cluster_taxon_lowest)
        plt.savefig(plot_out_mth_current_cluster, dpi=300)
        plt.close('all')

cluster_to_genome_member_file_handle.close()
cluster_to_taxon_file_handle.close()

