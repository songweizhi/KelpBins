import os


def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)
    return folder_list


Ribosomal_related_COG_file = '/Users/songweizhi/COG_DB/Ribosomal_protein_related_COGs.txt'
COG_annotation_wd_pcofg = '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_pcofg_COG_wd'


# get ribosomal protein related COGs
Ribosomal_related_COGs = set()
for each in open(Ribosomal_related_COG_file):
    each_split = each.strip().split(' ')
    Ribosomal_related_COGs.add(each_split[1])

#
all_faa_wd_list = get_no_hidden_folder_list(COG_annotation_wd_pcofg)


for each_folder in all_faa_wd_list:
    each_folder_base_name = '_'.join(each_folder.split('_')[:-2])
    protein_id_cog_file = '%s/%s/%s_protein-id_cog.txt' % (COG_annotation_wd_pcofg, each_folder, each_folder_base_name)
    non_ribosomal_cog_num = 0
    ribosomal_cog_num = 0
    for each_gene_cog in open(protein_id_cog_file):
        each_gene_cog_split = each_gene_cog.strip().split('\t')
        #print(each_gene_cog_split)

        if each_gene_cog_split[1] in Ribosomal_related_COGs:
            ribosomal_cog_num += 1
        else:
            non_ribosomal_cog_num += 1

    non_ribosomal_cog_percent = float("{0:.2f}".format(non_ribosomal_cog_num*100/(non_ribosomal_cog_num + ribosomal_cog_num)))
    ribosomal_cog_percent = float("{0:.2f}".format(ribosomal_cog_num*100/(non_ribosomal_cog_num + ribosomal_cog_num)))

    print('%s\t%s\t%s' % (each_folder_base_name, non_ribosomal_cog_percent, ribosomal_cog_percent))

