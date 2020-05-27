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


wd = '/Users/songweizhi/Desktop/KelpBins'
os.chdir(wd)

HGT_PG_validated_pcofg =          'combined_pcofg/PG_pcofg_new_right_algorithm/GoodBins_0.5_0.05_PG_pcofg_normal.txt'
PG_validated_pcofg_gene_id_cog =  '/Users/songweizhi/Desktop/000/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_protein-id_cog.txt'
PG_validated_pcofg_cog_stats =    '/Users/songweizhi/Desktop/000/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_cog_stats.txt'
PG_validated_pcofg_gene_id_kegg = '/Users/songweizhi/Desktop/000/GoodBins_0.5_0.05_PG_pcofg_recipient_gene_KO_assignment_DCBA.txt'
Flanking_region_plots_pcofg =     'combined_pcofg/Flanking_region_plots_pcofg'
recent_hgt_iden_cutoff = 95

recent_hgt_function_file =                'combined_pcofg/GoodBins_0.5_0.05_pcofg_recent_hgt_%s fun.txt' % recent_hgt_iden_cutoff
Flanking_region_plots_pcofg_recent_hgt =  'combined_pcofg/GoodBins_0.5_0.05_pcofg_recent_hgt_%s_flk_plots' % recent_hgt_iden_cutoff


# create folder
force_create_folder(Flanking_region_plots_pcofg_recent_hgt)


gene_to_ko_id_dict = {}
gene_to_ko_desc_dict = {}
for each_gene_KO in open(PG_validated_pcofg_gene_id_kegg):

    if not each_gene_KO.startswith('Gene_id'):
        each_gene_KO_split = each_gene_KO.strip().split('\t')
        gene_id = each_gene_KO_split[0]
        KO_id = 'NA'
        KO_desc = 'NA'
        if len(each_gene_KO_split) > 1:
            KO_id = each_gene_KO_split[1]
            KO_desc = each_gene_KO_split[5]
        gene_to_ko_id_dict[gene_id] = KO_id
        gene_to_ko_desc_dict[KO_id] = KO_desc


gene_to_cog_dict = {}
for each_gene_COG in open(PG_validated_pcofg_gene_id_cog):
    each_gene_COG_split = each_gene_COG.strip().split('\t')
    gene_id= each_gene_COG_split[0]
    COG_id= each_gene_COG_split[1]
    gene_to_cog_dict[gene_id] = COG_id


cog_id_to_description_dict = {}
for each_cog_id in open(PG_validated_pcofg_cog_stats):
    each_cog_id_split = each_cog_id.strip().split('\t')
    cog_id = each_cog_id_split[0]
    cog_description = each_cog_id_split[1]
    cog_id_to_description_dict[cog_id] = cog_description


recent_hgt_function_file_handle = open(recent_hgt_function_file, 'w')
num = 0
for hgt_o in open(HGT_PG_validated_pcofg):
    if not hgt_o.startswith('Gene_1'):
        hgt_o_split = hgt_o.strip().split('\t')
        gene_1 = hgt_o_split[0]
        gene_2 = hgt_o_split[1]

        concatenated_gene_1_2 = '%s___%s' % (gene_1, gene_2)

        if concatenated_gene_1_2 not in suspicious_HGTs:

            identity = float(hgt_o_split[2])
            gene_1_cog = 'NA'
            gene_2_cog = 'NA'
            gene_1_description = 'NA'
            gene_2_description = 'NA'

            gene_1_ko = 'NA'
            if gene_1 in gene_to_ko_id_dict:
                gene_1_ko = gene_to_ko_id_dict[gene_1]

            gene_2_ko = 'NA'
            if gene_2 in gene_to_ko_id_dict:
                gene_2_ko = gene_to_ko_id_dict[gene_2]


            gene_1_ko_desc = 'NA'
            if gene_1_ko in gene_to_ko_desc_dict:
                gene_1_ko_desc = gene_to_ko_desc_dict[gene_1_ko]

            gene_2_ko_desc = 'NA'
            if gene_2_ko in gene_to_ko_desc_dict:
                gene_2_ko_desc = gene_to_ko_desc_dict[gene_2_ko]

            if gene_1 in gene_to_cog_dict:
                gene_1_cog = gene_to_cog_dict[gene_1]
            if gene_2 in gene_to_cog_dict:
                gene_2_cog = gene_to_cog_dict[gene_2]
            if gene_1_cog in cog_id_to_description_dict:
                gene_1_description = cog_id_to_description_dict[gene_1_cog]
            if gene_2_cog in cog_id_to_description_dict:
                gene_2_description = cog_id_to_description_dict[gene_2_cog]

            for_out_gene_1 = '%s\t%s\t%s\t%s\t%s\n' % (gene_1, gene_1_cog, gene_1_description, gene_1_ko, gene_1_ko_desc)
            for_out_gene_2 = '%s\t%s\t%s\t%s\t%s\n' % (gene_2, gene_2_cog, gene_2_description, gene_2_ko, gene_2_ko_desc)

            if identity >= recent_hgt_iden_cutoff:
                recent_hgt_function_file_handle.write(for_out_gene_1)
                recent_hgt_function_file_handle.write(for_out_gene_2)
                num += 2

                # get flanking region plot for recent hgts
                flk_plot = '%s___%s.eps' % (gene_1, gene_2)
                pwd_flk_plot = '%s/%s' % (Flanking_region_plots_pcofg, flk_plot)
                os.system('cp %s %s/' % (pwd_flk_plot, Flanking_region_plots_pcofg_recent_hgt))


print(num)
recent_hgt_function_file_handle.close()