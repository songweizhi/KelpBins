import os
import shutil
from copy import deepcopy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


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


def extract_gene_flanking_sequences(prodigal_output_folder, bin_file_no_ext, flk_extraction_gene_list, flanking_length, min_flanking_length, output_seq_file_handle):

    pwd_gbk_file = '%s/%s.gbk' % (prodigal_output_folder, bin_file_no_ext)

    for contig in SeqIO.parse(pwd_gbk_file, 'genbank'):
        contig_id = contig.id
        contig_length = len(contig.seq)

        for gene in contig.features:

            if 'locus_tag' in gene.qualifiers:
                gene_id = gene.qualifiers['locus_tag'][0]

                if gene_id in flk_extraction_gene_list:
                    gene_start = gene.location.start
                    gene_end = gene.location.end
                    gene_left_flk_start = gene_start - flanking_length
                    gene_left_flk_end = gene_start
                    gene_right_flk_start = gene_end
                    gene_right_flk_end = gene_end + flanking_length


                    if gene_left_flk_start < 1:
                        gene_left_flk_start = 1

                    if gene_right_flk_end > contig_length:
                        gene_right_flk_end = contig_length

                    gene_left_flk_length = gene_left_flk_end - gene_left_flk_start
                    gene_right_flk_length = gene_right_flk_end - gene_right_flk_start

                    gene_left_flk_seq = str(contig.seq)[gene_left_flk_start - 1 : gene_left_flk_end - 1]
                    gene_right_flk_seq = str(contig.seq)[gene_right_flk_start:gene_right_flk_end]

                    # write out left flk seq
                    if gene_left_flk_length >= min_flanking_length:
                        output_seq_file_handle.write('>%s_L%s\n' % (gene_id, gene_left_flk_length))
                        output_seq_file_handle.write('%s\n' % (gene_left_flk_seq))

                    # write out right flk seq
                    if gene_right_flk_length >= min_flanking_length:
                        output_seq_file_handle.write('>%s_R%s\n' % (gene_id, gene_right_flk_length))
                        output_seq_file_handle.write('%s\n' % (gene_right_flk_seq))


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


bin_folder =                    '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05'
detected_HGTs =                 '/Users/songweizhi/Desktop/KelpBins/combined_pcofg/GoodBins_0.5_0.05_PG_pcofg_normal.txt'
cluster_to_genome_member_file = '/Users/songweizhi/Desktop/KelpBins/genome_culster_summary_folder/cluster_to_genome_member.txt'
prodigal_output_folder =        '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_all_prodigal_output'
multi_level_detection =         True
flanking_length =               1000
min_flanking_length =           30
distance_cutoff =               300


# output folder
HGT_with_time_wd = '/Users/songweizhi/Desktop/HGT_with_time_wd'


force_create_folder(HGT_with_time_wd)


# read in cluster info
cluster_to_genome_member_dict = {}
genome_to_cluster_dict = {}
n = 0
bin_num_list = []
for bin_cluster in open(cluster_to_genome_member_file):
    bin_cluster_split = bin_cluster.strip().split(':\t')
    cluster_ID = bin_cluster_split[0]
    bin_member = bin_cluster_split[1].split('\t')

    bin_num_list.append(len(bin_member))

    if len(bin_member) > 8:
        print(cluster_ID)
        print(bin_member)
        print(len(bin_member))
        print()
        n += 1

# print(n)
#
# print(bin_num_list)
# print(sorted(bin_num_list))
#
