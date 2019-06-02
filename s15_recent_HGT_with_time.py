import os
import shutil
from s0_suspicious_HGTs import suspicious_HGTs


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


recent_HGT_iden =           95
HGT_flk_folder =            '/Users/songweizhi/Desktop/KelpBins/combined_pcofg/Flanking_region_plots_pcofg_manually_checked'
recent_HGT_flk_folder =     '/Users/songweizhi/Desktop/KelpBins/combined_pcofg/GoodBins_0.5_0.05_pcofg_recent_hgt_95_flk_plots'


force_create_folder(recent_HGT_flk_folder)


recent_recipient_set = set()
for each_HGT in open('/Users/songweizhi/Desktop/GoodBins_0.5_0.05_PG_pcofg_normal.txt'):

    if not each_HGT.startswith('Gene_1'):
        each_HGT_split = each_HGT.strip().split('\t')
        gene_1 = each_HGT_split[0]
        gene_2 = each_HGT_split[1]
        identity = float(each_HGT_split[2])
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        direction = each_HGT_split[6]

        recipient_genome = direction.split('-->')[1]
        if '%)' in recipient_genome:
            recipient_genome = recipient_genome.split('(')[0]

        if gene_1_genome == recipient_genome:
            recipient_gene = gene_1
        else:
            recipient_gene = gene_2

        if (identity >= recent_HGT_iden) and ('%s___%s' % (gene_1, gene_2) not in suspicious_HGTs):
            recent_recipient_set.add(recipient_gene)
            pwd_flk_plot = '%s/%s___%s.eps' % (HGT_flk_folder, gene_1, gene_2)
            os.system('cp %s %s/' % (pwd_flk_plot, recent_HGT_flk_folder))

print(recent_recipient_set)


# get the cluster of recent recipients
genome_to_cluster_dict = {}
cluster_to_genome_number_dict = {}
for each_cluster in open('/Users/songweizhi/Desktop/KelpBins/genome_culster_summary_folder/cluster_to_genome_member.txt'):
    cluster_ID = each_cluster.strip().split(':\t')[0]
    cluster_ID_genome_member = each_cluster.strip().split(':\t')[1].split('\t')
    cluster_ID_genome_member_no_ext = ['.'.join(i.split('.')[:-1]) for i in cluster_ID_genome_member]
    for genome in cluster_ID_genome_member_no_ext:
        genome_to_cluster_dict[genome] = cluster_ID

    cluster_to_genome_number_dict[cluster_ID] = len(cluster_ID_genome_member_no_ext)


for recent_recipient in recent_recipient_set:

    recent_recipient_genome = '_'.join(recent_recipient.split('_')[:-1])
    recent_recipient_genome_cluster = genome_to_cluster_dict[recent_recipient_genome]

    for_print = '%s\t%s\t%s' % (recent_recipient, recent_recipient_genome_cluster, cluster_to_genome_number_dict[recent_recipient_genome_cluster])

    print(for_print)




