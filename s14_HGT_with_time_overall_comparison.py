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
pwd_combined_528bins_db =       '/Users/songweizhi/Desktop/KelpBins/combined_528bins/combined_528bins.fasta'

multi_level_detection =         True
flanking_length =               1000
min_flanking_length =           30
distance_cutoff =               50

# output folder
HGT_with_time_wd = '/Users/songweizhi/Desktop/HGT_with_time_wd'


force_create_folder(HGT_with_time_wd)


# read in cluster info
cluster_to_genome_member_dict = {}
genome_to_cluster_dict = {}
for bin_cluster in open(cluster_to_genome_member_file):
    bin_cluster_split = bin_cluster.strip().split(':\t')
    cluster_ID = bin_cluster_split[0]
    bin_member = bin_cluster_split[1].split('\t')

    # get cluster_to_genome_member_dict
    if len(bin_member) > 1:
        cluster_to_genome_member_dict[cluster_ID] = bin_member

    # get genome_to_cluster_dict
    for each_bin in bin_member:
        each_bin_no_ext = '.'.join(each_bin.split('.')[:-1])
        genome_to_cluster_dict[each_bin_no_ext] = cluster_ID


# get genome_to_HGT and cluster_to_HGT dict
recipient_genome_to_HGT_dict = {}
cluster_to_HGT_dict = {}
for each_HGT in open(detected_HGTs):
    if not each_HGT.startswith('Gene_1'):
        each_HGT_split = each_HGT.strip().split('\t')
        gene_1 = each_HGT_split[0]
        gene_2 = each_HGT_split[1]
        genome_1 = '_'.join(gene_1.split('_')[:-1])
        genome_2 = '_'.join(gene_2.split('_')[:-1])

        direction = each_HGT_split[5]
        if multi_level_detection is True:
            direction = each_HGT_split[6]

        if '%)' in direction:
            direction = direction.split('(')[0]

        if genome_1 == direction.split('-->')[1]:
            recipient_genome = genome_1
            recipient_gene = gene_1
        else:
            recipient_gene = gene_2
            recipient_genome = genome_2

        # get cluster of recipient genome
        recipient_genome_cluster = genome_to_cluster_dict[recipient_genome]

        # get recipient_genome_to_HGT_dict
        if recipient_genome not in recipient_genome_to_HGT_dict:
            recipient_genome_to_HGT_dict[recipient_genome] = [recipient_gene]
        else:
            recipient_genome_to_HGT_dict[recipient_genome].append(recipient_gene)

        # get cluster_to_HGT_dict
        if recipient_genome_cluster not in cluster_to_HGT_dict:
            cluster_to_HGT_dict[recipient_genome_cluster] = [recipient_gene]
        else:
            cluster_to_HGT_dict[recipient_genome_cluster].append(recipient_gene)


for each_cluster in cluster_to_genome_member_dict:
    each_cluster_genome_list = cluster_to_genome_member_dict[each_cluster]

    each_cluster_HGTs_num = 0
    if each_cluster in cluster_to_HGT_dict:
        each_cluster_HGTs_num = len(cluster_to_HGT_dict[each_cluster])

    if each_cluster_HGTs_num > 0:

        current_cluster_HGT_time_wd =           '%s/%s'                         % (HGT_with_time_wd, each_cluster)
        current_cluster_HGT_flk_seqs =          '%s/%s_flk_seqs.fasta'          % (current_cluster_HGT_time_wd, each_cluster)
        current_cluster_HGT_flk_blast =          '%s/%s_flk_blastn.tab'         % (current_cluster_HGT_time_wd, each_cluster)

        os.mkdir(current_cluster_HGT_time_wd)

        current_cluster_HGT_flk_seqs_handle = open(current_cluster_HGT_flk_seqs, 'w')
        for each_genome in each_cluster_genome_list:
            genome_name_no_ext = '.'.join(each_genome.split('.')[:-1])

            # combine current cluster genomes and add genome name to contig ID
            pwd_genome = '%s/%s' % (bin_folder, each_genome)

            # HGT from current genome
            current_genome_HGTs = []
            if genome_name_no_ext in recipient_genome_to_HGT_dict:
                current_genome_HGTs = recipient_genome_to_HGT_dict[genome_name_no_ext]

            #print('%s\t%s\t%s' % (each_cluster, each_genome, '\t'.join(current_genome_HGTs)))

            if current_genome_HGTs != []:
                extract_gene_flanking_sequences(prodigal_output_folder, genome_name_no_ext, current_genome_HGTs, flanking_length, min_flanking_length, current_cluster_HGT_flk_seqs_handle)

        current_cluster_HGT_flk_seqs_handle.close()


        # compare between extracted flanking sequences and genome sequences
        blastn_cmd = 'blastn -query %s -db %s -outfmt 6 -out %s' % (current_cluster_HGT_flk_seqs, pwd_combined_528bins_db, current_cluster_HGT_flk_blast)
        os.system(blastn_cmd)


        print(each_cluster)
        get_distance_between_matched_blocks(current_cluster_HGT_flk_blast, distance_cutoff)
        print()

