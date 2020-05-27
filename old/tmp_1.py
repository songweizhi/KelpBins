import os


def Get_circlize_plot(multi_level_detection, output_prefix, pwd_candidates_file_PG_normal_txt, genome_to_taxon_dict, circos_HGT_R, pwd_plot_circos, taxon_rank, taxon_rank_num, pwd_MetaCHIP_op_folder):

    pwd_cir_plot_t1 =              '%s/%s_%s_cir_plot_t1.txt'              % (pwd_MetaCHIP_op_folder, output_prefix, taxon_rank_num)
    pwd_cir_plot_t1_sorted =       '%s/%s_%s_cir_plot_t1_sorted.txt'       % (pwd_MetaCHIP_op_folder, output_prefix, taxon_rank_num)
    pwd_cir_plot_t1_sorted_count = '%s/%s_%s_cir_plot_t1_sorted_count.txt' % (pwd_MetaCHIP_op_folder, output_prefix, taxon_rank_num)
    pwd_cir_plot_matrix_filename = '%s/%s_%s_cir_plot_matrix.csv'          % (pwd_MetaCHIP_op_folder, output_prefix, taxon_rank_num)


    name2taxon_dict = {}
    transfers = []
    for each in open(pwd_candidates_file_PG_normal_txt):
        if not each.startswith('Gene_1'):
            each_split = each.strip().split('\t')
            Gene_1 = each_split[0]
            Gene_2 = each_split[1]
            Genome_1 = '_'.join(Gene_1.split('_')[:-1])
            Genome_2 = '_'.join(Gene_2.split('_')[:-1])

            if Genome_1 in genome_to_taxon_dict:
                Genome_1_taxon = genome_to_taxon_dict[Genome_1]
            else:
                Genome_1_taxon = '%s_' % taxon_rank


            if Genome_2 in genome_to_taxon_dict:
                Genome_2_taxon = genome_to_taxon_dict[Genome_2]
            else:
                Genome_2_taxon = '%s_' % taxon_rank

            Direction = each_split[5]
            if multi_level_detection == True:
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
    if len(all_group_id) == 1:
        print('Too less group (1), plot skipped')
    elif 1 < len(all_group_id) <= 200:
        os.system('Rscript %s -m %s -p %s' % (circos_HGT_R, pwd_cir_plot_matrix_filename, pwd_plot_circos))
    else:
        print('Too many groups (>200), plot skipped')

    # rm tmp files
    # os.system('rm %s' % pwd_cir_plot_t1)
    # os.system('rm %s' % pwd_cir_plot_t1_sorted)
    # os.system('rm %s' % pwd_cir_plot_t1_sorted_count)


multi_level_detection = True
output_prefix = 'Test'
taxon_rank = 'p'
taxon_rank_num = 10
wd =                                '/Users/songweizhi/Desktop/ring_chart'
pwd_plot_circos =                   '%s/Grouping_Figure_1.svg'                                              % wd
taxon_info =                 '%s/Kelp_NM_dRep99_p32_grouping.txt'                                              % wd

circos_HGT_R =                      '%s/circos_HGT_with_subgroup.R' % wd

pwd_candidates_file_PG_normal_txt = '%s/interested_hgts.txt' % wd
pwd_grouping_file =                        '%s/bin_bin_grouping.txt'              % wd




MAG_to_taxon_dict = {}
for genome in open(taxon_info):
    MAG_id = genome.strip().split(',')[1]
    taxon_name = genome.strip().split(',')[2][3:]
    MAG_to_taxon_dict[MAG_id] = taxon_name


taxon_to_group_id_dict = {}
for group in open(pwd_grouping_file):
    group_id = group.strip().split(',')[0]
    group_taxon = group.strip().split(',')[2]
    if group_id not in taxon_to_group_id_dict:
        taxon_to_group_id_dict[group_id] = group_taxon

genome_to_taxon_dict = {}
for genome in open(pwd_grouping_file):
    group_id2 = genome.strip().split(',')[0]
    genome_name = genome.strip().split(',')[1]
    genome_to_taxon_dict[genome_name] = taxon_to_group_id_dict[group_id2]


Get_circlize_plot(multi_level_detection,
                  output_prefix,
                  pwd_candidates_file_PG_normal_txt,
                  genome_to_taxon_dict,
                  circos_HGT_R,
                  pwd_plot_circos,
                  taxon_rank,
                  taxon_rank_num,
                  wd)


cir_plot_matrix =               '%s/Test_10_cir_plot_matrix.csv'                % wd
cir_plot_matrix_with_phylum =   '%s/Test_10_cir_plot_matrix_with_phylum.csv'    % wd


cir_plot_matrix_with_phylum_handle = open(cir_plot_matrix_with_phylum, 'w')
n = 0
for line in open(cir_plot_matrix):
    line_split = line.strip().split('\t')
    if n == 0:
        header_list = [ '%s__%s' % (MAG_to_taxon_dict[i], i) for i in line_split]
        cir_plot_matrix_with_phylum_handle.write('\t%s\n' % '\t'.join(header_list))
    else:
        cir_plot_matrix_with_phylum_handle.write('%s__%s' % (MAG_to_taxon_dict[line_split[0]], line))
    n += 1
cir_plot_matrix_with_phylum_handle.close()


os.system('Rscript %s/circos_HGT_with_subgroup.R -m %s/Test_10_cir_plot_matrix_with_phylum.csv -p %s/Test_10_cir_plot_matrix_with_phylum.pdf -s __' % (wd, wd, wd))


'''

cd /Users/songweizhi/Desktop/SpongeEMP/interested_HGTs/test
Rscript SpongeEMP_circos_HGT_with_subgroup.R -m Test_10_cir_plot_matrix_with_phylum.csv -p Test_10_cir_plot_matrix_with_phylum.pdf -s __

'''

