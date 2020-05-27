import os


def get_kelp_mag_taxon(grouping_file, grouping_kelp_only, identified_taxon):
    grouping_kelp_only_handle = open(grouping_kelp_only, 'w')
    identified_taxon_handle = open(identified_taxon, 'w')
    identified_taxon_list = []
    for mag in open(grouping_file):

        mag_split = mag.strip().split(',')

        if '_Refined_' in mag_split[1]:
            grouping_kelp_only_handle.write('%s\t%s\n' % (mag_split[1], mag_split[2]))

            if mag_split[2] not in identified_taxon_list:
                identified_taxon_handle.write('%s\n' % mag_split[2])
                identified_taxon_list.append(mag_split[2])

    grouping_kelp_only_handle.close()
    identified_taxon_handle.close()


# file in
grouping_file_folder = '/Users/songweizhi/Desktop/Kelp_NM/Figure_data/Figure_Tree_with_functions/Kelp/grouping_files'
grouping_p = '%s/Kelp_NM_dRep99_p32_grouping.txt'       % grouping_file_folder
grouping_c = '%s/Kelp_NM_dRep99_c50_grouping.txt'       % grouping_file_folder
grouping_o = '%s/Kelp_NM_dRep99_o129_grouping.txt'      % grouping_file_folder


# file out
grouping_p_kelp_only =      '%s/Kelp_dRep99_p_grouping.txt'         % grouping_file_folder
grouping_c_kelp_only =      '%s/Kelp_dRep99_c_grouping.txt'         % grouping_file_folder
grouping_o_kelp_only =      '%s/Kelp_dRep99_o_grouping.txt'         % grouping_file_folder
identified_taxon_p =        '%s/Kelp_dRep99_p_taxons.txt'           % grouping_file_folder
identified_taxon_c =        '%s/Kelp_dRep99_c_taxons.txt'           % grouping_file_folder
identified_taxon_o =        '%s/Kelp_dRep99_o_taxons.txt'           % grouping_file_folder
grouping_p_kelp_only_iTOL = '%s/Kelp_dRep99_p_grouping_iTOL.txt'    % grouping_file_folder
grouping_c_kelp_only_iTOL = '%s/Kelp_dRep99_c_grouping_iTOL.txt'    % grouping_file_folder
grouping_o_kelp_only_iTOL = '%s/Kelp_dRep99_o_grouping_iTOL.txt'    % grouping_file_folder


# get input files for BioSAK's iTOL module
get_kelp_mag_taxon(grouping_p, grouping_p_kelp_only, identified_taxon_p)
get_kelp_mag_taxon(grouping_c, grouping_c_kelp_only, identified_taxon_c)
get_kelp_mag_taxon(grouping_o, grouping_o_kelp_only, identified_taxon_o)


# run BioSAK
BioSAK_cmd_p = 'BioSAK iTOL -ColorRange -LeafGroup %s -LegendTitle Phylum -out %s' % (grouping_p_kelp_only, grouping_p_kelp_only_iTOL)
BioSAK_cmd_c = 'BioSAK iTOL -ColorStrip -LeafGroup %s -LegendTitle Class -out %s'  % (grouping_c_kelp_only, grouping_c_kelp_only_iTOL)
BioSAK_cmd_o = 'BioSAK iTOL -ColorStrip -LeafGroup %s -LegendTitle Order -out %s'  % (grouping_o_kelp_only, grouping_o_kelp_only_iTOL)

os.system(BioSAK_cmd_p)
os.system(BioSAK_cmd_c)
os.system(BioSAK_cmd_o)




