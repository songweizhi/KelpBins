import Kelp_HGT_config
from Kelp_HGT_config import Kelp_bin_list
from Kelp_HGT_config import Tara_NM_bin_list


# get Kelp and Tara MAG total size
Kelp_bin_total_size = 0
Tara_NM_bin_total_size = 0
for genome in open(Kelp_HGT_config.Kelp_NM_dRep99_genome_size_file):
    if not genome.startswith('Genome	Size(Mbp)'):
        genome_id = genome.strip().split('\t')[0].split('.fasta')[0]
        genome_size = float(genome.strip().split('\t')[1])
        if genome_id in Kelp_bin_list:
            Kelp_bin_total_size += genome_size
        if genome_id in Tara_NM_bin_list:
            Tara_NM_bin_total_size += genome_size


Kelp_bin_recipient_gene_list = []
Tara_NM_bin_recipient_gene_list = []
for each_HGT in open(Kelp_HGT_config.Kelp_NM_dRep99_pcofg_detected_HGTs):

    if not each_HGT.startswith('Gene_1'):
        each_HGT_split = each_HGT.strip().split('\t')

        each_HGT_gene_1 = each_HGT_split[0]
        each_HGT_gene_2 = each_HGT_split[1]
        each_HGT_genome_1 = '_'.join(each_HGT_gene_1.split('_')[:-1])
        each_HGT_genome_2 = '_'.join(each_HGT_gene_2.split('_')[:-1])
        concatenated_genes = '%s___%s' % (each_HGT_gene_1, each_HGT_gene_2)
        direction = each_HGT_split[6]

        if '%)' in direction:
            direction = direction.split('(')[0]

        if each_HGT_genome_1 == direction.split('-->')[1]:
            recipient_genome = each_HGT_genome_1
            recipient_gene = each_HGT_gene_1
        else:
            recipient_gene = each_HGT_gene_2
            recipient_genome = each_HGT_genome_2

        if recipient_genome in Kelp_bin_list:
            Kelp_bin_recipient_gene_list.append(recipient_gene)

        if recipient_genome in Tara_NM_bin_list:
            Tara_NM_bin_recipient_gene_list.append(recipient_gene)


######################################################## report ########################################################

print('\n================ Report ===================\n')

print('\t\t\t\t\tKelp\t\tTara_NM')
print('Bin number:\t\t\t%s\t\t\t %s' % (len(Kelp_bin_list), len(Tara_NM_bin_list)))
print('Total size (Mbp):\t%s\t\t %s' % (float("{0:.2f}".format(Kelp_bin_total_size)), float("{0:.2f}".format(Tara_NM_bin_total_size))))
print('HGT number:\t\t\t%s\t\t\t %s' % (len(Kelp_bin_recipient_gene_list), len(Tara_NM_bin_recipient_gene_list)))
print('HGT frequency:\t\t%s\t\t %s'  % ((float("{0:.2f}".format(len(Kelp_bin_recipient_gene_list)/Kelp_bin_total_size))), float("{0:.2f}".format(len(Tara_NM_bin_recipient_gene_list)/Tara_NM_bin_total_size))))

print('\n===========================================\n')

