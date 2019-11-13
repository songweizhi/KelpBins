import os


def Get_circlize_plot(multi_level_detection, output_prefix, pwd_candidates_file_PG_normal_txt, circos_HGT_R, pwd_plot_circos, pwd_MetaCHIP_op_folder):

    pwd_cir_plot_t1 =              '%s/%s_cir_plot_t1.txt'              % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_t1_sorted =       '%s/%s_cir_plot_t1_sorted.txt'       % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_t1_sorted_count = '%s/%s_cir_plot_t1_sorted_count.txt' % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_matrix_filename = '%s/%s_cir_plot_matrix.csv'          % (pwd_MetaCHIP_op_folder, output_prefix)

    transfers = []
    for each in open(pwd_candidates_file_PG_normal_txt):
        if not each.startswith('Gene_1'):
            each_split = each.strip().split('\t')

            Direction = each_split[5]
            if multi_level_detection == True:
                Direction = each_split[6]

            if '%)' in Direction:
                Direction = Direction.split('(')[0]

            transfers.append(Direction)

    tmp1 = open(pwd_cir_plot_t1, 'w')
    all_location_id = []
    for each_t in transfers:
        each_t_split = each_t.split('-->')
        donor = each_t_split[0]
        recipient = each_t_split[1]
        donor_location = donor.split('_')[0]
        recipient_location = recipient.split('_')[0]

        if donor_location not in all_location_id:
            all_location_id.append(donor_location)
        if recipient_location not in all_location_id:
            all_location_id.append(recipient_location)
        tmp1.write('%s,%s\n' % (donor_location, recipient_location))
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

    all_location_id = sorted(all_location_id)

    matrix_file = open(pwd_cir_plot_matrix_filename, 'w')
    matrix_file.write('\t' + '\t'.join(all_location_id) + '\n')
    for each_1 in all_location_id:
        row = [each_1]
        for each_2 in all_location_id:
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


multi_level_detection =             True
output_prefix =                     'Test'
pwd_candidates_file_PG_normal_txt = '/Users/songweizhi/Desktop/KelpBins/between_location/GoodBins_0.5_0.05_PG_pcofg_normal.txt'
pwd_plot_circos =                   '/Users/songweizhi/Desktop/KelpBins/between_location/GoodBins_0.5_0.05_PG_pcofg_normal_location.png'
circos_HGT_R =                      '/Users/songweizhi/Desktop/KelpBins/between_location/MetaCHIP_circos_HGT.R'
pwd_MetaCHIP_op_folder =            '/Users/songweizhi/Desktop/KelpBins/between_location'


Get_circlize_plot(multi_level_detection, output_prefix, pwd_candidates_file_PG_normal_txt, circos_HGT_R, pwd_plot_circos, pwd_MetaCHIP_op_folder)

