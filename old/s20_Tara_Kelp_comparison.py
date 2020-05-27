
genome_size_file =              '/Users/songweizhi/Desktop/Tara_Kelp_all_genome_size.txt'
pwd_detected_HGT_txt_combined = '/Users/songweizhi/Desktop/Tara_Kelp_pcofg_detected_HGTs.txt'


Kelp_genome_list = []
Tara_genome_list = []
genome_size_dict = {}
Kelp_genome_total_size = 0
Tara_genome_total_size = 0
for each_genome_size in open(genome_size_file):
    if not each_genome_size.startswith('Genome'):
        each_genome_size_split = each_genome_size.strip().split('\t')
        genome_id = each_genome_size_split[0]
        genome_size = float(each_genome_size_split[1])
        genome_size_dict[genome_id] = genome_size

        if '_Refined_' in genome_id:
            Kelp_genome_list.append(genome_id)
            Kelp_genome_total_size += genome_size
        else:
            Tara_genome_list.append(genome_id)
            Tara_genome_total_size += genome_size


Kelp_genome_HGT_num = 0
Tara_genome_HGT_num = 0
for each in open(pwd_detected_HGT_txt_combined):
    if not each.startswith('Gene_1'):

        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        direction = each_split[6]
        plot_file = '%s___%s.SVG' % (gene_1, gene_2)

        recipient_genome = direction.split('-->')[1]
        if '%)' in recipient_genome:
            recipient_genome = recipient_genome.split('(')[0]

        if '_Refined_' in recipient_genome:
            Kelp_genome_HGT_num += 1
        else:
            Tara_genome_HGT_num += 1


Kelp_genome_mean_size = float("{0:.2f}".format(Kelp_genome_total_size/len(Kelp_genome_list)))
Tara_genome_mean_size = float("{0:.2f}".format(Tara_genome_total_size/len(Tara_genome_list)))


print('Kelp genome mean size: %s Mbp' % Kelp_genome_mean_size)
print('Tara genome mean size: %s Mbp' % Tara_genome_mean_size)
print('Kelp genome HGT num: %s'       % Kelp_genome_HGT_num)
print('Tara genome HGT num: %s'       % Tara_genome_HGT_num)
print('Kelp genome HGT/Mbp %s'        % float("{0:.2f}".format(Kelp_genome_HGT_num/Kelp_genome_total_size)))
print('Tara genome HGT/Mbp %s'        % float("{0:.2f}".format(Tara_genome_HGT_num/Tara_genome_total_size)))

