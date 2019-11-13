import os
import glob
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
import s0_Kelp_bins_config
import seaborn as sns
from s0_suspicious_HGTs import suspicious_HGTs


def Get_circlize_plot(output_prefix, pwd_candidates_file_PG_normal_txt, genome_to_taxon_dict, circos_HGT_R, pwd_plot_circos, taxon_rank):

    pwd_cir_plot_t1 =              '%s_%s_cir_plot_t1.txt'              % (output_prefix, taxon_rank_num)
    pwd_cir_plot_t1_sorted =       '%s_%s_cir_plot_t1_sorted.txt'       % (output_prefix, taxon_rank_num)
    pwd_cir_plot_t1_sorted_count = '%s_%s_cir_plot_t1_sorted_count.txt' % (output_prefix, taxon_rank_num)
    pwd_cir_plot_matrix_filename = '%s_%s_cir_plot_matrix.csv'          % (output_prefix, taxon_rank_num)


    name2taxon_dict = {}
    transfers = []
    for each in open(pwd_candidates_file_PG_normal_txt):
        if not each.startswith('Gene_1'):
            each_split = each.strip().split('\t')
            Gene_1 = each_split[0]
            Gene_2 = each_split[1]
            Genome_1 = '_'.join(Gene_1.split('_')[:-1])
            Genome_2 = '_'.join(Gene_2.split('_')[:-1])
            #Genome_1_taxon = genome_to_taxon_dict[Genome_1]
            #Genome_2_taxon = genome_to_taxon_dict[Genome_2]

            if Genome_1 in genome_to_taxon_dict:
                Genome_1_taxon = genome_to_taxon_dict[Genome_1]
            else:
                Genome_1_taxon = '%s_' % taxon_rank


            if Genome_2 in genome_to_taxon_dict:
                Genome_2_taxon = genome_to_taxon_dict[Genome_2]
            else:
                Genome_2_taxon = '%s_' % taxon_rank


            Direction = each_split[6]
            if '%)' in Direction:
                Direction = Direction.split('(')[0]

            if Genome_1 not in name2taxon_dict:
                name2taxon_dict[Genome_1] = Genome_1_taxon
            if Genome_2 not in name2taxon_dict:
                name2taxon_dict[Genome_2] = Genome_2_taxon
            transfers.append(Direction)


    tmp1 = open(pwd_cir_plot_t1, 'w')
    all_group_id = []
    for each_t in transfers:
        each_t_split = each_t.split('-->')
        donor = each_t_split[0]
        recipient = each_t_split[1]
        donor_id = name2taxon_dict[donor]
        recipient_id = name2taxon_dict[recipient]
        if donor_id not in all_group_id:
            all_group_id.append(donor_id)
        if recipient_id not in all_group_id:
            all_group_id.append(recipient_id)
        tmp1.write('%s,%s\n' % (donor_id, recipient_id))
    tmp1.close()

    os.system('cat %s | sort > %s' % (pwd_cir_plot_t1, pwd_cir_plot_t1_sorted))

    current_t = ''
    count = 0
    tmp2 = open(pwd_cir_plot_t1_sorted_count, 'w')
    for each_t2 in open(pwd_cir_plot_t1_sorted):
        each_t2 = each_t2.strip()
        if current_t == '':
            current_t = each_t2
            count += 1
        elif current_t == each_t2:
            count += 1
        elif current_t != each_t2:
            tmp2.write('%s,%s\n' % (current_t, count))
            current_t = each_t2
            count = 1
    tmp2.write('%s,%s\n' % (current_t, count))
    tmp2.close()

    # read in count as dict
    transfer_count = {}
    for each_3 in open(pwd_cir_plot_t1_sorted_count):
        each_3_split = each_3.strip().split(',')
        key = '%s,%s' % (each_3_split[0], each_3_split[1])
        value = each_3_split[2]
        transfer_count[key] = value

    all_group_id = sorted(all_group_id)

    matrix_file = open(pwd_cir_plot_matrix_filename, 'w')
    matrix_file.write('\t' + '\t'.join(all_group_id) + '\n')
    for each_1 in all_group_id:
        row = [each_1]
        for each_2 in all_group_id:
            current_key = '%s,%s' % (each_2, each_1)
            if current_key not in transfer_count:
                row.append('0')
            else:
                row.append(transfer_count[current_key])
        matrix_file.write('\t'.join(row) + '\n')
    matrix_file.close()

    # get plot with R
    os.system('Rscript %s -m %s -p %s' % (circos_HGT_R, pwd_cir_plot_matrix_filename, pwd_plot_circos))

    # rm tmp files
    os.system('rm %s' % pwd_cir_plot_t1)
    os.system('rm %s' % pwd_cir_plot_t1_sorted)
    os.system('rm %s' % pwd_cir_plot_t1_sorted_count)



####################################################### file I/O #######################################################

wd = '/Users/songweizhi/Desktop/KelpBins/combined_pcofg'
grouping_file_folder =                  'GoodBins_0.5_0.05_pcofg_grouping'
PG_output_folder =                      'TT_90MGs_PG'
pwd_candidates_seq_file =               'GoodBins_0.5_0.05_all_combined_ffn.fasta'
output_prefix =                         'GoodBins_0.5_0.05'
detection_ranks =                       'pcofg'
grouping_file_highest_rank =            '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_pcofg_grouping/GoodBins_0.5_0.05_p10_grouping.txt'
pwd_candidates_file_PG_normal_txt =     '%s_PG_%s_normal.txt'   % (output_prefix, detection_ranks)
plot_all_genomes =                      False


########################################################################################################################

os.chdir(wd)
circos_HGT_R =                          '/Users/songweizhi/PycharmProjects/MetaCHIP/MetaCHIP/MetaCHIP_circos_HGT.R'


grouping_file_re = '%s/%s_*_grouping.txt' % (s0_Kelp_bins_config.pwd_grouping_file_folder, output_prefix)
PG_grouping_file_list_with_path = [('%s/%s' % (s0_Kelp_bins_config.pwd_grouping_file_folder, os.path.basename(file_name))) for file_name in glob.glob(grouping_file_re)]

# get input genome list
input_genome_list = set()
for each_genome in open(grouping_file_highest_rank):
    input_genome_list.add(each_genome.strip().split(',')[1])


# store genome size in dict
genome_size_dict = {}
for each_size in open(s0_Kelp_bins_config.genome_size_file):
    if each_size.strip() != 'Genome\tSize(Mbp)':
        each_size_split = each_size.strip().split('\t')
        genome_file_name = each_size_split[0]
        genome_size_Mbp = float(each_size_split[1])
        genome_size_dict[genome_file_name] = genome_size_Mbp


# get recipient gene list
recipient_genome_list = []
for each in open(pwd_candidates_file_PG_normal_txt):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        concatenated_genes = '%s___%s' % (gene_1, gene_2)

        if concatenated_genes not in suspicious_HGTs:
            gene_1_genome = '_'.join(gene_1.split('_')[:-1])
            gene_2_genome = '_'.join(gene_2.split('_')[:-1])
            direction = each_split[6]
            recipient_genome = direction.split('-->')[1]
            if '%)' in recipient_genome:
                recipient_genome = recipient_genome.split('(')[0]
            recipient_genome_list.append(recipient_genome)


# get the number of HGT detected from each genome
genome_HGT_num_absolute_dict = {}
for input_genome in input_genome_list:
    genome_HGT_num_absolute_dict[input_genome] = recipient_genome_list.count(input_genome)


genome_HGT_num_norm_dict = {}
for genome in genome_HGT_num_absolute_dict:
    genome_size = genome_size_dict[genome]
    HGT_num_per_Mbp = float("{0:.2f}".format(genome_HGT_num_absolute_dict[genome]/genome_size))
    genome_HGT_num_norm_dict[genome] = HGT_num_per_Mbp


for grouping_file in PG_grouping_file_list_with_path:

    grouping_file_name = grouping_file.strip().split('/')[-1]
    taxon_rank_num = grouping_file_name[len(output_prefix)+1:].split('_')[0]
    taxon_rank = taxon_rank_num[0]

    pwd_group_to_taxon_file = '%s/%s_%s_group_to_taxon.txt' % (s0_Kelp_bins_config.pwd_grouping_file_folder, output_prefix, taxon_rank_num)

    pwd_plot = '%s_%s_HGT_num_plot.png' % (output_prefix, taxon_rank_num)
    pwd_plot_circos = '%s_%s_HGT_circos.png' % (output_prefix, taxon_rank_num)

    taxon_to_group_id_dict = {}
    for group in open(pwd_group_to_taxon_file):
        group_id = group.strip().split(',')[0]
        group_taxon = group.strip().split(',')[1]
        taxon_to_group_id_dict[group_id] = group_taxon


    # get genome to taxon dict
    genome_to_taxon_dict = {}
    for genome in open(grouping_file):
        group_id2 = genome.strip().split(',')[0]
        genome_name = genome.strip().split(',')[1]
        genome_to_taxon_dict[genome_name] = taxon_to_group_id_dict[group_id2]


    genome_HGT_num_label_lol_all = []
    genome_HGT_num_label_lol_non_zero = []
    for genome2 in input_genome_list:

        HGT_num_absolute = genome_HGT_num_absolute_dict[genome2]
        HGT_num_norm = genome_HGT_num_norm_dict[genome2]

        if genome2 in genome_to_taxon_dict:
            genome_taxon = genome_to_taxon_dict[genome2]
        else:
            genome_taxon = '%s_' % taxon_rank

        if HGT_num_absolute != 0:
            genome_HGT_num_label_lol_non_zero.append([genome2, HGT_num_absolute, HGT_num_norm, genome_taxon])
        genome_HGT_num_label_lol_all.append([genome2, HGT_num_absolute, HGT_num_norm, genome_taxon])


    genome_HGT_num_label_lol_all_df = pd.DataFrame(np.array(genome_HGT_num_label_lol_all, dtype=object), columns=['genome', 'absolute', 'norm', 'label'])

    genome_HGT_num_label_lol_non_zero_df = pd.DataFrame(np.array(genome_HGT_num_label_lol_non_zero, dtype=object), columns=['genome', 'absolute', 'norm', 'label'])

    if plot_all_genomes == True:
        df_to_plot = genome_HGT_num_label_lol_all_df
    else:
        df_to_plot = genome_HGT_num_label_lol_non_zero_df


    # get plot
    sns.set(style="whitegrid")  # whitegrid
    ax = sns.boxplot(x="label", y="norm", data=df_to_plot, orient='v', width=0.5, linewidth=0.5, showfliers=False)
    ax = sns.stripplot(x="label", y="norm", data=df_to_plot, color='orange', size=2, orient='v', jitter=0.23)
    ax.set(xlabel='Taxon', ylabel='Number of HGT/Mbp sequence')

    plt.xticks(rotation=270)

    # save plot
    plt.tight_layout()
    plt.savefig(pwd_plot, dpi=300)
    plt.clf()

    Get_circlize_plot(output_prefix, pwd_candidates_file_PG_normal_txt, genome_to_taxon_dict, circos_HGT_R, pwd_plot_circos, taxon_rank)
