import Kelp_HGT_config


def get_HGT_frequency(detected_HGT_file, genome_size_file):

    # get genome total size
    genome_total_size = 0
    genome_num = 0
    for genome in open(genome_size_file):
        if not genome.startswith('Genome	Size(Mbp)'):
            genome_size = float(genome.strip().split('\t')[1])
            genome_total_size += genome_size
            genome_num += 1

    # get total_number of detected HGTs
    detected_HGT_num = len(open(detected_HGT_file).readlines()) - 1

    # get HGT frequency
    HGT_frequency = detected_HGT_num / genome_total_size
    HGT_frequency = (float("{0:.2f}".format(HGT_frequency)))

    return detected_HGT_num, genome_num, (float("{0:.2f}".format(genome_total_size))), HGT_frequency


Kelp_all_detected_HGT_num, Kelp_all_genome_num, Kelp_all_genome_total_size, Kelp_all_HGT_frequency =        get_HGT_frequency(Kelp_HGT_config.Kelp_bins_detected_HGTs, Kelp_HGT_config.genome_size_Kelp_bins)
Kelp_dRep_detected_HGT_num, Kelp_dRep_genome_num, Kelp_dRep_genome_total_size, Kelp_dRep_HGT_frequency =    get_HGT_frequency(Kelp_HGT_config.Kelp_dereplicated_bins_detected_HGTs, Kelp_HGT_config.genome_size_Kelp_dereplicated_bins)
Tara_NM_detected_HGT_num, Tara_NM_genome_num, Tara_NM_genome_total_size, Tara_NM_HGT_frequency =            get_HGT_frequency(Kelp_HGT_config.Tara_NM_detected_HGTs, Kelp_HGT_config.genome_size_Tara_NM_bins)
Tara_SD_detected_HGT_num, Tara_SD_genome_num, Tara_SD_genome_total_size, Tara_SD_HGT_frequency =            get_HGT_frequency(Kelp_HGT_config.Tara_SD_detected_HGTs, Kelp_HGT_config.genome_size_Tara_SD_bins)


print('Dataset\tBinNum\tSize\tHgtNum\tFrequency')
print('Kelp_all\t%s\t%s\t%s\t%s'  % (Kelp_all_genome_num,  Kelp_all_genome_total_size,  Kelp_all_detected_HGT_num,  Kelp_all_HGT_frequency))
print('Kelp_dRep\t%s\t%s\t%s\t%s' % (Kelp_dRep_genome_num, Kelp_dRep_genome_total_size, Kelp_dRep_detected_HGT_num, Kelp_dRep_HGT_frequency))
print('Tara_NM\t%s\t%s\t%s\t%s'   % (Tara_NM_genome_num,   Tara_NM_genome_total_size,   Tara_NM_detected_HGT_num,   Tara_NM_HGT_frequency))
print('Tara_SD\t%s\t%s\t%s\t%s'   % (Tara_SD_genome_num,   Tara_SD_genome_total_size,   Tara_SD_detected_HGT_num,   Tara_SD_HGT_frequency))


HGT_all_bin_list = set()
for HGT_all_bin in open(Kelp_HGT_config.Kelp_bins_detected_HGTs):
    if not HGT_all_bin.startswith('Gene_1'):
        HGT_all_bin_split = HGT_all_bin.strip().split('\t')
        gene_1 = HGT_all_bin_split[0]
        gene_2 = HGT_all_bin_split[1]
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        direction = HGT_all_bin_split[6]
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

        HGT_all_bin_list.add(recipient_gene)


HGT_dRep_bin_list = set()
for HGT_all_bin in open(Kelp_HGT_config.Kelp_dereplicated_bins_detected_HGTs):
    if not HGT_all_bin.startswith('Gene_1'):
        HGT_all_bin_split = HGT_all_bin.strip().split('\t')
        gene_1 = HGT_all_bin_split[0]
        gene_2 = HGT_all_bin_split[1]
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        direction = HGT_all_bin_split[6]
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

        HGT_dRep_bin_list.add(recipient_gene)


print(HGT_all_bin_list)
print(len(HGT_all_bin_list))

print(HGT_dRep_bin_list)
print(len(HGT_dRep_bin_list))


m = 0
n = 0
for HGT in HGT_dRep_bin_list:

    if HGT not in HGT_all_bin_list:
        n += 1
    else:
        m += 1

print()
print(m)
print(n)




