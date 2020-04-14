

def get_taxon_stats(GTDB_output_file, summary_rank):

    taxon_assignment_dict = {}
    for each_genome in open(GTDB_output_file):
        if not each_genome.startswith('user_genome'):
            each_split = each_genome.strip().split('\t')
            bin_name = each_split[0]

            assignment_full = []
            if len(each_split) == 1:
                assignment_full = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
            elif (len(each_split) > 1) and (';' in each_split[1]):
                assignment = each_split[1].split(';')
                if len(assignment) == 7:
                    assignment_full = assignment
                if len(assignment) == 6:
                    assignment_full = assignment + ['s__']
                if len(assignment) == 5:
                    assignment_full = assignment + ['g__', 's__']
                if len(assignment) == 4:
                    assignment_full = assignment + ['f__', 'g__', 's__']
                if len(assignment) == 3:
                    assignment_full = assignment + ['o__', 'f__', 'g__', 's__']
                if len(assignment) == 2:
                    assignment_full = assignment + ['c__', 'o__', 'f__', 'g__', 's__']

            elif (len(each_split) > 1) and (';' not in each_split[1]):
                assignment_full = [each_split[1]] + ['p__', 'c__', 'o__', 'f__', 'g__', 's__']

            # store in dict
            taxon_assignment_dict[bin_name] = assignment_full

    # get all identified taxon at defined ranks
    rank_to_position_dict = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
    specified_rank_pos = rank_to_position_dict[summary_rank]
    identified_taxon_list = []
    for each_TaxonAssign in taxon_assignment_dict:
        specified_rank_id = taxon_assignment_dict[each_TaxonAssign][specified_rank_pos]
        if specified_rank_id not in identified_taxon_list:
            identified_taxon_list.append(specified_rank_id)

    taxon_2_genome_dict = {}
    for each_taxon in identified_taxon_list:
        genome_list = []
        for genome in taxon_assignment_dict:
            if taxon_assignment_dict[genome][specified_rank_pos] == each_taxon:
                genome_list.append(genome)
        taxon_2_genome_dict[each_taxon] = genome_list

    taxon_2_genome_num_dict = {}
    for i in taxon_2_genome_dict:
        taxon_2_genome_num_dict[i] = len(taxon_2_genome_dict[i])

    return taxon_2_genome_num_dict


wd = '/Users/songweizhi/Desktop/Kelp_NM/Taxon_comparison'

Kelp_GTDB_txt = '%s/Kelp_GTDB_r89.tsv'      % wd
Tara_GTDB_txt = '%s/Tara_NM_GTDB_r89.tsv'   % wd
summary_rank =  'p'

Kelp_taxon_stats_dict = get_taxon_stats(Kelp_GTDB_txt, summary_rank)
Tara_taxon_stats_dict = get_taxon_stats(Tara_GTDB_txt, summary_rank)


identified_taxon_list = []
for i in Kelp_taxon_stats_dict:
    identified_taxon_list.append(i)

for i in Tara_taxon_stats_dict:
    if i not in identified_taxon_list:
        identified_taxon_list.append(i)

for i in sorted(identified_taxon_list):

    i_num_in_kelp = 0
    if i in Kelp_taxon_stats_dict:
        i_num_in_kelp = Kelp_taxon_stats_dict[i]

    i_num_in_Tara = 0
    if i in Tara_taxon_stats_dict:
        i_num_in_Tara = Tara_taxon_stats_dict[i]

    print('%s\t%s\t%s' % (i, i_num_in_kelp, i_num_in_Tara))


