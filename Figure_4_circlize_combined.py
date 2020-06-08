import os
import glob
import itertools
from Bio import SeqIO
from string import ascii_uppercase


def get_group_index_list():

    def iter_all_strings():
        size = 1
        while True:
            for s in itertools.product(ascii_uppercase, repeat=size):
                yield "".join(s)
            size += 1

    group_index_list = []
    for s in iter_all_strings():
        group_index_list.append(s)
        if s == 'ZZZ':
            break

    return group_index_list


def Get_circlize_plot(multi_level_detection, output_prefix, pwd_candidates_file_PG_normal_txt, genome_to_taxon_dict, circos_HGT_R, pwd_plot_circos, taxon_rank, pwd_MetaCHIP_op_folder):

    pwd_cir_plot_t1 =              '%s/%s_cir_plot_t1.txt'              % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_t1_sorted =       '%s/%s_cir_plot_t1_sorted.txt'       % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_t1_sorted_count = '%s/%s_cir_plot_t1_sorted_count.txt' % (pwd_MetaCHIP_op_folder, output_prefix)
    pwd_cir_plot_matrix_filename = '%s/%s_cir_plot_matrix.csv'          % (pwd_MetaCHIP_op_folder, output_prefix)


    name2taxon_dict = {}
    transfers = []
    for each in open(pwd_candidates_file_PG_normal_txt):
        if not each.startswith('Gene_1'):
            each_split      = each.strip().split('\t')
            Gene_1          = each_split[0]
            Gene_2          = each_split[1]
            Genome_1        = '_'.join(Gene_1.split('_')[:-1])
            Genome_2        = '_'.join(Gene_2.split('_')[:-1])
            Genome_1_taxon  = genome_to_taxon_dict[Genome_1]
            Genome_2_taxon  = genome_to_taxon_dict[Genome_2]

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


# file in
wd                                  = '/Users/songweizhi/Desktop/ring_chart/combined'
detected_kelp_HGTs_seq              = '/Users/songweizhi/Desktop/Kelp_NM/0_file_in/Kelp_NM_dRep99_pcofg_detected_HGTs_recipient_genes.faa'
detected_HGTs                       = '/Users/songweizhi/Desktop/Kelp_NM/0_file_in/Kelp_NM_dRep99_pcofg_detected_HGTs.txt'
gene_annotation_folder              = '/Users/songweizhi/Desktop/Kelp_NM/Figure_data/Figure_Tree_with_functions/Kelp_and_Tara/COG_KEGG_dbCAN/COG_KEGG_dbCAN_combined'
interested_functions                = ['COG1653', 'COG1175', 'K02027', 'COG0601', 'K02035', 'COG0410', 'K01996', 'CE6', 'COG3119', 'COG0346', 'K02610', 'K02611', 'K02613', 'K01912', 'K02618', 'COG4663', 'K03320', 'COG0004', 'COG2128', 'COG0376', 'K03782', 'COG2128', 'COG0714', 'K01626', 'K01995', 'K02031', 'K02032', 'K02033', 'K01996', 'K01997', 'K01999', 'K02034', 'K03076', 'K02035', 'K01897', 'K10914', 'K01998', 'COG0051', 'COG0100', 'COG0048']
transporters                        = ['COG1653', 'COG1175', 'K02027', 'COG0601', 'K02035', 'COG0410', 'K01996']
sugar_and_phlorotannins_degradation = ['CE6', 'COG3119', 'COG0346', 'K02610', 'K02611', 'K02613', 'K01912', 'K02618', 'COG4663', 'K03320', 'COG0004']
ROS_and_stress_response             = ['COG2128', 'COG0376', 'K03782', 'COG2128', 'COG0714']
quorum_sensing                      = ['K01626', 'K01995', 'K02031', 'K02032', 'K02033', 'K01996', 'K01997', 'K01999', 'K02034', 'K03076', 'K02035', 'K01897', 'K10914', 'K01998']
ribosomal_proteins                  = ['COG0051', 'COG0100', 'COG0048']
output_prefix                       = 'combined'
taxon_rank                          = 'o'
multi_level_detection               = True

pwd_grouping_file = '%s/bin_bin_grouping.txt'                       % wd
if taxon_rank ==        'p':
    taxon_info_raw =    '%s/Kelp_NM_dRep99_p32_grouping.txt'        % '/Users/songweizhi/Desktop/ring_chart'
    taxon_info =        '%s/Kelp_NM_dRep99_p32_grouping.txt'        % wd
    circos_HGT_R =      '%s/circos_HGT_with_subgroup_p.R'           % wd
if taxon_rank ==        'c':
    taxon_info_raw =    '%s/Kelp_NM_dRep99_c50_grouping.txt'        % '/Users/songweizhi/Desktop/ring_chart'
    taxon_info =        '%s/Kelp_NM_dRep99_c50_grouping.txt'        % wd
    circos_HGT_R =      '%s/circos_HGT_with_subgroup_c.R'           % wd
if taxon_rank ==        'o':
    taxon_info_raw =    '%s/Kelp_NM_dRep99_o129_grouping.txt'       % '/Users/songweizhi/Desktop/ring_chart'
    taxon_info =        '%s/Kelp_NM_dRep99_o129_grouping.txt'       % wd
    circos_HGT_R =      '%s/circos_HGT_with_subgroup_o.R'           % wd


# file out
interested_hgts =                       '%s/interested_hgts.txt'                        % wd
cir_plot_matrix =                       '%s/%s_cir_plot_matrix.csv'                     % (wd, output_prefix)
cir_plot_matrix_with_rank =             '%s/%s_cir_plot_matrix_%s.csv'                  % (wd, output_prefix, taxon_rank)
pwd_plot_circos =                       '%s/%s_Grouping_Figure_1.svg'                   % (wd, output_prefix)
bin_bin_grouping_file =                 '%s/bin_bin_grouping.txt'                       % wd


gene_annotation_folder_re = '%s/*.txt' % gene_annotation_folder
gene_annotation_file_list = [os.path.basename(file_name) for file_name in glob.glob(gene_annotation_folder_re)]

gene_to_function_dict = {}
for gene_annotation_file in gene_annotation_file_list:
    pwd_gene_annotation_file = '%s/%s' % (gene_annotation_folder, gene_annotation_file)
    if '_Refined_' in gene_annotation_file:
        for gene in open(pwd_gene_annotation_file):
            gene_split = gene.strip().split('\t')
            gene_to_function_dict[gene_split[0]] = gene_split[1:]

# get function id to cate dict
function_id_to_cate_dict = {}
for i in transporters:
    if i not in function_id_to_cate_dict:
        function_id_to_cate_dict[i] = ['transporters']
for i in sugar_and_phlorotannins_degradation:
    if i not in function_id_to_cate_dict:
        function_id_to_cate_dict[i] = ['sugar_and_phlorotannins_degradation']
    else:
        function_id_to_cate_dict[i].append('sugar_and_phlorotannins_degradation')
for i in ROS_and_stress_response:
    if i not in function_id_to_cate_dict:
        function_id_to_cate_dict[i] = ['ROS_and_stress_response']
    else:
        function_id_to_cate_dict[i].append('ROS_and_stress_response')
for i in quorum_sensing:
    if i not in function_id_to_cate_dict:
        function_id_to_cate_dict[i] = ['quorum_sensing']
    else:
        function_id_to_cate_dict[i].append('quorum_sensing')
for i in ribosomal_proteins:
    if i not in function_id_to_cate_dict:
        function_id_to_cate_dict[i] = ['ribosomal_proteins']
    else:
        function_id_to_cate_dict[i].append('ribosomal_proteins')

interested_recipient_list = []

transporters_num = 0
sugar_and_phlorotannins_degradation_num = 0
ROS_and_stress_response_num = 0
quorum_sensing_num = 0
ribosomal_proteins_num = 0
for seq_record in SeqIO.parse(detected_kelp_HGTs_seq, 'fasta'):
    if '_Refined_' in seq_record.id:
        seq_record_functions = gene_to_function_dict[seq_record.id]
        interested_hgt = 'No'
        for gene_function in seq_record_functions:
            if gene_function in interested_functions:
                interested_hgt = 'Yes'

            if gene_function in transporters:
                transporters_num += 1
            if gene_function in sugar_and_phlorotannins_degradation:
                sugar_and_phlorotannins_degradation_num += 1
            if gene_function in ROS_and_stress_response:
                ROS_and_stress_response_num += 1
            if gene_function in quorum_sensing:
                quorum_sensing_num += 1
            if gene_function in ribosomal_proteins:
                ribosomal_proteins_num += 1

        if interested_hgt == 'Yes':
            interested_recipient_list.append(seq_record.id)


print(transporters_num)
print(sugar_and_phlorotannins_degradation_num)
print(ROS_and_stress_response_num)
print(quorum_sensing_num)
print(ribosomal_proteins_num)
print(len(interested_recipient_list))

interested_hgts_handle = open(interested_hgts, 'w')
genome_to_fun_cate_dict = {}
all_genome_set = set()
n = 0
for each in open(detected_HGTs):
    if not each.startswith('Gene_1'):

        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        direction = each_split[6]

        recipient_genome = direction.split('-->')[1]
        if '%)' in recipient_genome:
            recipient_genome = recipient_genome.split('(')[0]

        donor_gene = ''
        donor_genome = ''
        recipient_gene = ''
        if gene_1_genome == recipient_genome:
            recipient_gene = gene_1
            donor_gene = gene_2
            donor_genome = gene_2_genome
        if gene_2_genome == recipient_genome:
            recipient_gene = gene_2
            donor_gene = gene_1
            donor_genome = gene_1_genome

        if ('_Refined_' in recipient_genome) and ('_Refined_' in donor_genome):

            if recipient_gene in interested_recipient_list:

                variant_genome = ''
                variant_genome_num = 0

                n += 1
                recipient_gene_function = gene_to_function_dict[recipient_gene]

                recipient_gene_function_cate_list = []
                for fun_id in recipient_gene_function:
                    if fun_id in function_id_to_cate_dict:
                        fun_id_cate = function_id_to_cate_dict[fun_id]
                        for individual_fun_id_cate in fun_id_cate:
                            if individual_fun_id_cate not in recipient_gene_function_cate_list:
                                recipient_gene_function_cate_list.append(individual_fun_id_cate)


                if donor_genome not in genome_to_fun_cate_dict:
                    genome_to_fun_cate_dict[donor_genome] = [recipient_gene_function_cate_list]

                    # write out without add index to the donor genome
                    interested_hgts_handle.write(each)
                    all_genome_set.add(donor_genome)
                    all_genome_set.add(recipient_genome)

                else:
                    if genome_to_fun_cate_dict[donor_genome] == [recipient_gene_function_cate_list]:

                        # write out without add index to the donor genome
                        interested_hgts_handle.write(each)
                        all_genome_set.add(donor_genome)
                        all_genome_set.add(recipient_genome)

                    else:

                        # check whether variant genome with the same function already exist
                        for gnm in genome_to_fun_cate_dict:

                            if (len(gnm.split('_')) == 6) and ('_'.join(gnm.split('_')[:5]) == donor_genome):

                                variant_genome_num += 1

                                if genome_to_fun_cate_dict[gnm] == [recipient_gene_function_cate_list]:
                                    variant_genome = gnm

                                    # rename to gnm and write out here
                                    if gene_1 == donor_gene:
                                        gene_1_variant = '%s_%s' % (gnm, gene_1.split('_')[-1])
                                        # interested_hgts_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s-->%s\n' % (gene_1_variant, gene_2, each_split[2], each_split[3], each_split[4], each_split[5], gnm, recipient_genome))
                                        interested_hgts_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s-->%s\t%s\n' % (gene_1_variant, gene_2, each_split[2], each_split[3], each_split[4], each_split[5], gnm, recipient_genome, '|'.join(recipient_gene_function_cate_list)))
                                    if gene_2 == donor_gene:
                                        gene_2_variant = '%s_%s' % (gnm, gene_2.split('_')[-1])
                                        interested_hgts_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s-->%s\t%s\n' % (gene_1, gene_2_variant, each_split[2], each_split[3], each_split[4], each_split[5], gnm, recipient_genome, '|'.join(recipient_gene_function_cate_list)))

                                    all_genome_set.add(gnm)
                                    all_genome_set.add(recipient_genome)

                        if variant_genome == '':
                            variant_genome = '%s_%s' % (donor_genome, variant_genome_num + 1)
                            genome_to_fun_cate_dict[variant_genome] = [recipient_gene_function_cate_list]
                            variant_genome_num += 1

                            # rename to variant_genome and write out here
                            if gene_1 == donor_gene:
                                gene_1_variant = '%s_%s' % (variant_genome, gene_1.split('_')[-1])
                                interested_hgts_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s-->%s\t%s\n' % (gene_1_variant, gene_2, each_split[2], each_split[3], each_split[4], each_split[5], variant_genome, recipient_genome, '|'.join(recipient_gene_function_cate_list)))
                            if gene_2 == donor_gene:
                                gene_2_variant = '%s_%s' % (variant_genome, gene_2.split('_')[-1])
                                interested_hgts_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s-->%s\t%s\n' % (gene_1, gene_2_variant, each_split[2], each_split[3], each_split[4], each_split[5], variant_genome, recipient_genome, '|'.join(recipient_gene_function_cate_list)))

                            all_genome_set.add(variant_genome)
                            all_genome_set.add(recipient_genome)


# get bin_bin_grouping file for all genomes
all_genome_list_sorted = sorted([i for i in all_genome_set])
group_index_list = get_group_index_list()
bin_bin_grouping_file_handle = open(bin_bin_grouping_file, 'w')
index = 0
for i in all_genome_list_sorted:
    bin_bin_grouping_file_handle.write('%s,%s,%s\n' % (group_index_list[index], i, i))
    index += 1
bin_bin_grouping_file_handle.close()


# update taxon file
taxon_dict = {}
for i in open(taxon_info_raw):
    i_split = i.strip().split(',')
    taxon_dict[i_split[1]] = [i_split[0], i_split[2]]

taxon_info_handle = open(taxon_info, 'w')
for i in all_genome_list_sorted:
    i_split = i.strip().split('_')
    if len(i_split) == 5:
        taxon_info_handle.write('%s,%s,%s\n' % (taxon_dict[i][0], i, taxon_dict[i][1]))
    if len(i_split) == 6:
        i_raw = '_'.join(i_split[:5])
        taxon_info_handle.write('%s,%s,%s\n' % (taxon_dict[i_raw][0], i, taxon_dict[i_raw][1]))
taxon_info_handle.close()


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


Get_circlize_plot(multi_level_detection, output_prefix, interested_hgts, genome_to_taxon_dict, circos_HGT_R, pwd_plot_circos, taxon_rank, wd)


cir_plot_matrix_with_phylum_handle = open(cir_plot_matrix_with_rank, 'w')
n = 0
for line in open(cir_plot_matrix):
    line_split = line.strip().split('\t')
    if n == 0:
        header_list = ['%s__%s' % (MAG_to_taxon_dict[i], i) for i in line_split]
        cir_plot_matrix_with_phylum_handle.write('\t%s\n' % '\t'.join(header_list))
    else:
        cir_plot_matrix_with_phylum_handle.write('%s__%s' % (MAG_to_taxon_dict[line_split[0]], line))
    n += 1
cir_plot_matrix_with_phylum_handle.close()


os.system('Rscript %s -m %s/%s_cir_plot_matrix_%s.csv -p %s/%s_cir_plot_matrix_%s.pdf -s __' % (circos_HGT_R, wd, output_prefix, taxon_rank, wd, output_prefix, taxon_rank))


############################################ set MAG color ############################################

fun_cate_to_color_dict = {'transporters':                                       '#DC7633',
                          'sugar_and_phlorotannins_degradation':                '#5DADE2',
                          'ROS_and_stress_response':                            '#2ECC71',
                          'quorum_sensing':                                     '#F1C40F',
                          'ribosomal_proteins':                                 '#A569BD',
                          'transporters__quorum_sensing':                       '#808000',
                          'ROS_and_stress_response__ROS_and_stress_response':   '#2ECC71'}

genome_to_color_dict = {}
for each in open(interested_hgts):
    if not each.startswith('Gene_1'):

        each_split = each.strip().split('\t')
        gene_1 = each_split[0]
        gene_2 = each_split[1]
        gene_1_genome = '_'.join(gene_1.split('_')[:-1])
        gene_2_genome = '_'.join(gene_2.split('_')[:-1])
        direction = each_split[6]

        recipient_genome = direction.split('-->')[1]
        if '%)' in recipient_genome:
            recipient_genome = recipient_genome.split('(')[0]

        donor_gene = ''
        donor_genome = ''
        recipient_gene = ''
        if gene_1_genome == recipient_genome:
            recipient_gene = gene_1
            donor_gene = gene_2
            donor_genome = gene_2_genome
        if gene_2_genome == recipient_genome:
            recipient_gene = gene_2
            donor_gene = gene_1
            donor_genome = gene_1_genome

        donor_gene_raw = donor_gene
        if len(donor_gene.split('_')) == 7:
            donor_gene_raw = '%s_%s' % ('_'.join(donor_gene.split('_')[:5]), donor_gene.split('_')[6])
        recipient_gene_raw = recipient_gene
        if len(recipient_gene.split('_')) == 7:
            recipient_gene_raw = '%s_%s' % ('_'.join(recipient_gene.split('_')[:5]), recipient_gene.split('_')[6])

        recipient_gene_fun = gene_to_function_dict[recipient_gene_raw]

        recipient_gene_fun_cate = []
        for i in recipient_gene_fun:
            if i in function_id_to_cate_dict:
                if recipient_gene_fun_cate == []:
                    recipient_gene_fun_cate = function_id_to_cate_dict[i]

        recipient_gene_fun_cate_color = ''
        if len(recipient_gene_fun_cate) == 1:
            recipient_gene_fun_cate_color = fun_cate_to_color_dict[recipient_gene_fun_cate[0]]
        if recipient_gene_fun_cate == ['ROS_and_stress_response', 'ROS_and_stress_response']:
            recipient_gene_fun_cate_color = '#3399FF'
        if recipient_gene_fun_cate == ['transporters', 'quorum_sensing']:
            recipient_gene_fun_cate_color = '#99CCFF'



        if recipient_genome not in genome_to_color_dict:
            genome_to_color_dict[recipient_genome] = ['grey90']
        else:
            if 'grey90' not in genome_to_color_dict[recipient_genome]:
                genome_to_color_dict[recipient_genome].append('grey90')

        if donor_genome not in genome_to_color_dict:
            genome_to_color_dict[donor_genome] = [recipient_gene_fun_cate_color]
        else:
            if recipient_gene_fun_cate_color not in genome_to_color_dict[donor_genome]:
                genome_to_color_dict[donor_genome].append(recipient_gene_fun_cate_color)

genome_to_color_dict['BH_ER_161216_Refined_10']   = ['grey90', '#99CCFF']
genome_to_color_dict['CB_ER_080217_Refined_25_2'] = ['#99CCFF']


for i in genome_to_color_dict:

    i_color_list = genome_to_color_dict[i]

    if len(i_color_list) == 1:
        print('%s__%s = \'%s\',' % (MAG_to_taxon_dict[i], i, i_color_list[0]))
    if len(i_color_list) == 2:
        if i_color_list[0] == 'grey90':
            print('%s__%s = \'%s\', # grey90' % (MAG_to_taxon_dict[i], i, i_color_list[1]))
        if i_color_list[1] == 'grey90':
            print('%s__%s = \'%s\', # grey90' % (MAG_to_taxon_dict[i], i, i_color_list[0]))
