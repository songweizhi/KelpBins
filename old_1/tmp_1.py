
def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def get_GTDB_LSA(query_genome_list, combined_GTDB_summary_files, full_lineage):

    # get genome_to_taxon_dict
    genome_to_taxon_dict = {}
    for genome_taxon in open(combined_GTDB_summary_files):
        if not genome_taxon.startswith('user_genome'):
            genome_taxon_split = genome_taxon.strip().split('\t')
            genome_to_taxon_dict[genome_taxon_split[0]] = genome_taxon_split[1]

    # get last common ancestor
    taxon_rank_list = ['d', 'p', 'c', 'o', 'f', 'g', 's']

    LSA_determined = False
    last_common_ancestor = None
    for taxon_rank in taxon_rank_list[::-1]:
        current_rank_index = taxon_rank_list.index(taxon_rank)

        current_rank_member_list = []
        for query_genome in query_genome_list:
            query_genome_taxon = genome_to_taxon_dict[query_genome].split(';')
            current_rank_member_list.append(query_genome_taxon[current_rank_index])

        current_rank_member_list_uniq = unique_list_elements(current_rank_member_list)

        if (LSA_determined == False) and (len(current_rank_member_list_uniq) == 1) and (current_rank_member_list_uniq[0] != '%s__' % taxon_rank):
            last_common_ancestor = current_rank_member_list_uniq[0]
            LSA_determined = True

        elif (full_lineage == True) and (LSA_determined == True):
            last_common_ancestor = current_rank_member_list_uniq[0] + ';' + last_common_ancestor

        if (taxon_rank == 'd') and (len(current_rank_member_list_uniq) > 1):
            last_common_ancestor = 'life'

    return last_common_ancestor



GTDB_results = '/Users/songweizhi/Desktop/KelpBins/gtdbtk.bac120.classification_528_r86.tsv'

query_genome_list = ['CB_ER_080217_Refined_38']




print(get_GTDB_LSA(query_genome_list, GTDB_results, full_lineage=True))


