import os
import glob


traits_to_plot = '''

pilus assmebly and adherence
COG4961, COG4965, COG4964 , K02651, K02279, K02280, K02282, K12510, K12511, K02278, K02283

conductance mechanosensitive channels (CMC)
COG1970 , K03282, K03442

aquaporins
K06188

transcriptional regulators
COG3682, COG1695, COG1510, COG3655, COG3279, K16137, K05800, K10778

signal transducers
COG0515, COG2972

fucoidans degradation
CAZy CBM47, GH29, GH95

alginate degradation
PL6, PL7, PL14, PL15,PL17

xylan degradation
CE3, CE6


A: Pilus assmebly and adherence

B: Conductance mechanosensitive channels
C: Aquaporins

D: Transcriptional regulators
E: Signal transducers

F: Fucoidans degradation
G: Alginate degradation
H: Xylan degradation

'''

pilus_assmebly_and_adherence_function_list  = ['COG4961', 'COG4965', 'COG4964', 'K02651', 'K02279', 'K02280', 'K02282', 'K12510', 'K12511', 'K02278', 'K02283']
CMC_function_list                           = ['COG1970', 'K03282', 'K03442']
aquaporins_function_list                    = ['K06188']
transcriptional_regulator_function_list     = ['COG3682', 'COG1695', 'COG1510', 'COG3655', 'COG3279', 'K16137', 'K05800', 'K10778']
signal_transducer_function_list             = ['COG0515', 'COG2972']
fucoidans_degradation_function_list         = ['CBM47', 'GH29', 'GH95']
alginate_degradation_function_list          = ['PL6', 'PL7', 'PL14', 'PL15', 'PL17']
xylan_degradation_function_list             = ['CE3', 'CE6']


# file in
annotation_reaults_folder                   = '/Users/songweizhi/Desktop/Kelp_NM/Figure_data/Figure_Tree_with_functions/COG_KEGG_dbCAN'
specified_function_list                     = xylan_degradation_function_list


file_out = ''
LegendTitle = ''
file_out_iTOL_compatible = ''
if specified_function_list == pilus_assmebly_and_adherence_function_list:
    file_out                    = '%s/pilus_assmebly_and_adherence_gene_proportion.txt' % annotation_reaults_folder
    LegendTitle                 = 'pilus_assmebly_and_adherence'
    file_out_iTOL_compatible    = '%s/pilus_assmebly_and_adherence_gene_proportion_iTOL.txt' % annotation_reaults_folder

if specified_function_list == CMC_function_list:
    file_out                    = '%s/CMC_gene_proportion.txt' % annotation_reaults_folder
    LegendTitle                 = 'CMC'
    file_out_iTOL_compatible    = '%s/CMC_gene_proportion_iTOL.txt' % annotation_reaults_folder

if specified_function_list == aquaporins_function_list:
    file_out                    = '%s/aquaporins_gene_proportion.txt' % annotation_reaults_folder
    LegendTitle                 = 'aquaporins'
    file_out_iTOL_compatible    = '%s/aquaporins_gene_proportion_iTOL.txt' % annotation_reaults_folder

if specified_function_list == transcriptional_regulator_function_list:
    file_out                    = '%s/transcriptional_regulator_gene_proportion.txt' % annotation_reaults_folder
    LegendTitle                 = 'transcriptional_regulator'
    file_out_iTOL_compatible    = '%s/transcriptional_regulator_gene_proportion_iTOL.txt' % annotation_reaults_folder

if specified_function_list == signal_transducer_function_list:
    file_out                    = '%s/signal_transducer_gene_proportion.txt' % annotation_reaults_folder
    LegendTitle                 = 'signal_transducer'
    file_out_iTOL_compatible    = '%s/signal_transducer_gene_proportion_iTOL.txt' % annotation_reaults_folder

if specified_function_list == fucoidans_degradation_function_list:
    file_out                    = '%s/fucoidans_degradation_gene_proportion.txt' % annotation_reaults_folder
    LegendTitle                 = 'fucoidans_degradation'
    file_out_iTOL_compatible    = '%s/fucoidans_degradation_gene_proportion_iTOL.txt' % annotation_reaults_folder

if specified_function_list == alginate_degradation_function_list:
    file_out                    = '%s/alginate_degradation_gene_proportion.txt' % annotation_reaults_folder
    LegendTitle                 = 'alginate_degradation'
    file_out_iTOL_compatible    = '%s/alginate_degradation_gene_proportion_iTOL.txt' % annotation_reaults_folder

if specified_function_list == xylan_degradation_function_list:
    file_out                    = '%s/xylan_degradation_gene_proportion.txt' % annotation_reaults_folder
    LegendTitle                 = 'xylan_degradation'
    file_out_iTOL_compatible    = '%s/xylan_degradation_gene_proportion_iTOL.txt' % annotation_reaults_folder


COG_annotation_results_re = '%s/COG/*.txt' % annotation_reaults_folder
COG_annotation_results_file_list = [os.path.basename(file_name) for file_name in glob.glob(COG_annotation_results_re)]
MAG_id_list = [i.split('_query_to_cog.txt')[0] for i in COG_annotation_results_file_list]


# combine annotations together
# for MAG in MAG_id_list:
#
#     pwd_MAG_COG_annotation_result   = '%s/COG/%s_query_to_cog.txt'     % (annotation_reaults_folder, MAG)
#     pwd_MAG_KEGG_annotation_result  = '%s/KEGG/%s_KO_assignment_D.txt' % (annotation_reaults_folder, MAG)
#     pwd_MAG_dbCAN_annotation_result = '%s/dbCAN/%s_dbCAN.txt'          % (annotation_reaults_folder, MAG)
#
#     current_MAG_gene_list = []
#
#     # get current_MAG_gene_to_COG_dict
#     current_MAG_gene_to_COG_dict = {}
#     for gene_cog in open(pwd_MAG_COG_annotation_result):
#         if not gene_cog.startswith('Query	COG	Category	Description'):
#             gene_cog_split = gene_cog.strip().split('\t')
#             gene_id = gene_cog_split[0]
#
#             if gene_id not in current_MAG_gene_list:
#                 current_MAG_gene_list.append(gene_id)
#
#             if len(gene_cog_split) > 1:
#                 cog_id = gene_cog_split[1]
#                 current_MAG_gene_to_COG_dict[gene_id] = cog_id
#
#     # get current_MAG_gene_to_KEGG_dict
#     current_MAG_gene_to_KEGG_dict = {}
#     for gene_kegg in open(pwd_MAG_KEGG_annotation_result):
#         gene_kegg_split = gene_kegg.strip().split('\t')
#         if len(gene_kegg_split) == 2:
#             gene_id = gene_kegg_split[0]
#             ko_id = gene_kegg_split[1]
#             current_MAG_gene_to_KEGG_dict[gene_id] = ko_id
#
#     # get current_MAG_gene_to_dbCAN_dict
#     current_MAG_gene_to_dbCAN_dict = {}
#     for gene_dbcan in open(pwd_MAG_dbCAN_annotation_result):
#         if not gene_dbcan.startswith('Query	Family	Activities'):
#             gene_dbcan_split = gene_dbcan.strip().split('\t')
#             gene_id = gene_dbcan_split[0]
#             cazy_id = gene_dbcan_split[1].split('.')[0]
#             current_MAG_gene_to_dbCAN_dict[gene_id] = cazy_id
#
#     # combine annotations together
#     pwd_COG_KEGG_dbCAN_combined_file = '%s/COG_KEGG_dbCAN_combined/%s_COG_KEGG_dbCAN.txt' % (annotation_reaults_folder, MAG)
#     pwd_COG_KEGG_dbCAN_combined_file_handle = open(pwd_COG_KEGG_dbCAN_combined_file, 'w')
#     for gene in current_MAG_gene_list:
#
#         current_gene_cog = 'noCOG'
#         if gene in current_MAG_gene_to_COG_dict:
#             current_gene_cog = current_MAG_gene_to_COG_dict[gene]
#
#         current_gene_kegg = 'noKEGG'
#         if gene in current_MAG_gene_to_KEGG_dict:
#             current_gene_kegg = current_MAG_gene_to_KEGG_dict[gene]
#
#         current_gene_CAZy = 'noCAZy'
#         if gene in current_MAG_gene_to_dbCAN_dict:
#             current_gene_CAZy = current_MAG_gene_to_dbCAN_dict[gene]
#
#         pwd_COG_KEGG_dbCAN_combined_file_handle.write('%s\t%s\t%s\t%s\n' % (gene, current_gene_cog, current_gene_kegg, current_gene_CAZy))
#
#     pwd_COG_KEGG_dbCAN_combined_file_handle.close()

file_out_handle = open(file_out, 'w')
for MAG in MAG_id_list:

    pwd_combined_annotation = '%s/COG_KEGG_dbCAN_combined/%s_COG_KEGG_dbCAN.txt' % (annotation_reaults_folder, MAG)

    gene_num_total = 0
    gene_num_specified_function = 0
    for each_gene in open(pwd_combined_annotation):

        each_gene_split = each_gene.strip().split('\t')
        each_gene_assigned_function_list = each_gene_split[1:]

        belong_to_specified_function = False
        for assigned_function in each_gene_assigned_function_list:
            if assigned_function in specified_function_list:
                belong_to_specified_function = True

        gene_num_total += 1
        if belong_to_specified_function is True:
            gene_num_specified_function += 1

    specified_function_gene_proportion = float("{0:.3f}".format(gene_num_specified_function*100/gene_num_total))
    file_out_handle.write('%s\t%s\n' % (MAG, specified_function_gene_proportion))

file_out_handle.close()

# get iTOL compatible file
BioSAK_cmd = 'BioSAK iTOL -Heatmap -LeafMatrix %s -LegendTitle %s -out %s' % (file_out, LegendTitle, file_out_iTOL_compatible)
os.system(BioSAK_cmd)

# BioSAK iTOL -Heatmap -LeafMatrix combined.txt -LegendTitle Proportion -out combined_iTOL.txt
