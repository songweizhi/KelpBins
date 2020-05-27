
def check_matched_element_in_list(to_check, list_in):

    matched = False
    for each in list_in:
        if to_check in each:
            matched = True

    return matched


def get_distance_between_matched_blocks(blastn_results, distance_cutoff):

    subject_to_query_dict = {}
    subject_to_query_gene_dict = {}
    best_blast_hit_dict = {}
    for blast_hit in open(blastn_results):

        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        query_gene = '_'.join(query.split('_')[:-1])
        query_genome = '_'.join(query_gene.split('_')[:-1])
        query_len = int(query.split('_')[-1][1:])
        subject = blast_hit_split[1]
        subject_genome = subject.split('___')[0]
        identity = float(blast_hit_split[2])
        align_len = int(blast_hit_split[3])
        query_cov = align_len/query_len
        bit_score = float(blast_hit_split[11])
        subject_start = min([int(blast_hit_split[8]), int(blast_hit_split[9])])
        subject_end = max([int(blast_hit_split[8]), int(blast_hit_split[9])])

        if (query_genome != subject_genome) and (query_cov >= 0.8):

            if subject not in subject_to_query_dict:
                subject_to_query_dict[subject] = [query]
                subject_to_query_gene_dict[subject] = {query_gene}
            else:
                subject_to_query_dict[subject].append(query)
                subject_to_query_gene_dict[subject].add(query_gene)

            info_to_keep = '%s_%s_____%s_____%s_____%s' % (query_gene, query.split('_')[-1][0], subject, subject_start, subject_end)
            match_key = '%s_%s_____%s' % (query_gene, query.split('_')[-1][0], subject)
            match_value = '%s_____%s_____%s' % (subject_start, subject_end, bit_score)

            if match_key not in best_blast_hit_dict:
                best_blast_hit_dict[match_key] = match_value
            else:
                if bit_score > float(best_blast_hit_dict[match_key].split('_____')[-1]):
                    best_blast_hit_dict[match_key] = match_value


    for each_subject in subject_to_query_gene_dict:

        current_subject_matched_gene = subject_to_query_gene_dict[each_subject]

        for gene in current_subject_matched_gene:

            gene_flk_l_re = '%s_L' % gene
            gene_flk_r_re = '%s_R' % gene

            if (check_matched_element_in_list(gene_flk_l_re, subject_to_query_dict[each_subject]) is True) and (check_matched_element_in_list(gene_flk_r_re, subject_to_query_dict[each_subject]) is True):

                left_side_match_key =  '%s_L_____%s' % (gene, each_subject)
                right_side_match_key = '%s_R_____%s' % (gene, each_subject)

                left_side_match_value = 'NA'
                if left_side_match_key in best_blast_hit_dict:
                    left_side_match_value = '-'.join(best_blast_hit_dict[left_side_match_key].split('_____')[:-1])

                right_side_match_value = 'NA'
                if right_side_match_key in best_blast_hit_dict:
                    right_side_match_value = '-'.join(best_blast_hit_dict[right_side_match_key].split('_____')[:-1])

                matched_blocks = [left_side_match_value, right_side_match_value]

                if int(right_side_match_value.split('-')[0]) > int(left_side_match_value.split('-')[1]):
                    matched_blocks_left_to_right = matched_blocks

                else:
                    matched_blocks_left_to_right = matched_blocks[::-1]

                distance_between_matched_blocks = int(matched_blocks_left_to_right[1].split('-')[0]) - int(matched_blocks_left_to_right[0].split('-')[1])

                if distance_between_matched_blocks < distance_cutoff:
                    print('%s\t%s\t%s\t%s\t%s' % (gene, each_subject, matched_blocks_left_to_right[0], matched_blocks_left_to_right[1], distance_between_matched_blocks))




blastn_results = '/Users/songweizhi/Desktop/HGT_with_time_wd/Cluster213/Cluster213_flk_blastn.tab'
distance_cutoff = 500




get_distance_between_matched_blocks(blastn_results, distance_cutoff)