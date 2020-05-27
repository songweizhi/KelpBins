import os
import glob
from Bio import SeqIO


# file in
detected_kelp_HGTs_seq              = '/Users/songweizhi/Desktop/Kelp_NM/0_file_in/Kelp_NM_dRep99_pcofg_detected_HGTs_recipient_genes.faa'
detected_HGTs                       = '/Users/songweizhi/Desktop/Kelp_NM/0_file_in/Kelp_NM_dRep99_pcofg_detected_HGTs.txt'
gene_annotation_folder              = '/Users/songweizhi/Desktop/Kelp_NM/Figure_data/Figure_Tree_with_functions/Kelp_and_Tara/COG_KEGG_dbCAN/COG_KEGG_dbCAN_combined'
interested_functions                = ['COG1653', 'COG1175', 'K02027', 'COG0601', 'K02035', 'COG0410', 'K01996', 'CE6', 'COG3119', 'COG0346', 'K02610', 'K02611', 'K02613', 'K01912', 'K02618', 'COG4663', 'K03320', 'COG0004', 'COG2128', 'COG0376', 'K03782', 'COG2128', 'COG0714', 'K01626', 'K01995', 'K02031', 'K02032', 'K02033', 'K01996', 'K01997', 'K01999', 'K02034', 'K03076', 'K02035', 'K01897', 'K10914', 'K01998', 'COG0051', 'COG0100', 'COG0048']
transporters                        = ['COG1653', 'COG1175', 'K02027', 'COG0601', 'K02035', 'COG0410', 'K01996']
sugar_and_phlorotannins_degradation = ['CE6', 'COG3119', 'COG0346', 'K02610', 'K02611', 'K02613', 'K01912', 'K02618', 'COG4663', 'K03320', 'COG0004']
ROS_and_stress_response             = ['COG2128', 'COG0376', 'K03782', 'COG2128', 'COG0714']
potentially_QS                      = ['K01626', 'K01995', 'K02031', 'K02032', 'K02033', 'K01996', 'K01997', 'K01999', 'K02034', 'K03076', 'K02035', 'K01897', 'K10914', 'K01998']
ribosomal_proteins                  = ['COG0051', 'COG0100', 'COG0048']

interested_functions = interested_functions

# file out
interested_hgts = '/Users/songweizhi/Desktop/ring_chart/interested_hgts.txt'


gene_annotation_folder_re = '%s/*.txt' % gene_annotation_folder
gene_annotation_file_list = [os.path.basename(file_name) for file_name in glob.glob(gene_annotation_folder_re)]


gene_to_function_dict = {}
for gene_annotation_file in gene_annotation_file_list:
    pwd_gene_annotation_file = '%s/%s' % (gene_annotation_folder, gene_annotation_file)
    if '_Refined_' in gene_annotation_file:
        for gene in open(pwd_gene_annotation_file):
            gene_split = gene.strip().split('\t')
            gene_to_function_dict[gene_split[0]] = gene_split[1:]


interested_recipient_list = []
for seq_record in SeqIO.parse(detected_kelp_HGTs_seq, 'fasta'):
    if '_Refined_' in seq_record.id:
        seq_record_functions = gene_to_function_dict[seq_record.id]
        interested_hgt = 'No'
        for gene_function in seq_record_functions:
            if gene_function in interested_functions:
                interested_hgt = 'Yes'
        if interested_hgt == 'Yes':
            interested_recipient_list.append(seq_record.id)


interested_hgts_handle = open(interested_hgts, 'w')
recipient_and_donor_genes = set()
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
                interested_hgts_handle.write(each)
                recipient_and_donor_genes.add(donor_gene)
                recipient_and_donor_genes.add(recipient_gene)


interested_hgts_handle.close()











id_to_fun_dict = {}

for fun_1 in transporters:
    id_to_fun_dict[fun_1] = ['transporters']

for fun_2 in sugar_and_phlorotannins_degradation:
    if fun_2 not in id_to_fun_dict:
        id_to_fun_dict[fun_2] = ['sugar_and_phlorotannins_degradation']
    else:
        id_to_fun_dict[fun_2].append('sugar_and_phlorotannins_degradation')

for fun_3 in ROS_and_stress_response:
    if fun_3 not in id_to_fun_dict:
        id_to_fun_dict[fun_3] = ['ROS_and_stress_response']
    else:
        id_to_fun_dict[fun_3].append('ROS_and_stress_response')

for fun_4 in potentially_QS:
    if fun_4 not in id_to_fun_dict:
        id_to_fun_dict[fun_4] = ['potentially_QS']
    else:
        id_to_fun_dict[fun_4].append('potentially_QS')

for fun_5 in ribosomal_proteins:
    if fun_5 not in id_to_fun_dict:
        id_to_fun_dict[fun_5] = ['ribosomal_proteins']
    else:
        id_to_fun_dict[fun_5].append('ribosomal_proteins')

print(recipient_and_donor_genes)
genome_to_function_dict= {}
for each_gene in recipient_and_donor_genes:

    each_gene_genome = '_'.join(each_gene.split('_')[:-1])
    each_gene_function = gene_to_function_dict[each_gene]

    each_gene_function_cate = set()

    for function in each_gene_function:
        if function in id_to_fun_dict:

            for i in id_to_fun_dict[function]:
                each_gene_function_cate.add(i)


    print(each_gene)
    print(each_gene_function)
    print(each_gene_function_cate)

    if each_gene_genome not in genome_to_function_dict:
        genome_to_function_dict[each_gene_genome] = each_gene_function_cate
    else:
        for i in each_gene_function_cate:
            genome_to_function_dict[each_gene_genome].add(i)

print(genome_to_function_dict)



for each in genome_to_function_dict:

    print(len(genome_to_function_dict[each]))







traits_to_plot = '''

1) transport
'COG1653', 'COG1175', 'K02027', 'COG0601', 'K02035', 'COG0410', 'K01996'

2) sugar and phlorotannins degradation
'CE6', 'COG3119', 'COG0346', 'K02610', 'K02611', 'K02613', 'K01912', 'K02618', 'COG4663', 'K03320', 'COG0004'

3) ROS and stress response
'COG2128', 'COG0376', 'K03782', 'COG2128', 'COG0714' 

4) potentially QS
'K01626', 'K01995', 'K02031', 'K02032', 'K02033', 'K01996', 'K01997', 'K01999', 'K02034', 'K03076', 'K02035', 'K01897', 'K10914', 'K01998'

5) ribosomal proteins
'COG0051', 'COG0100', 'COG0048'

'''
