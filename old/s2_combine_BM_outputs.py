import os
import glob
from Bio import SeqIO


def combine_BM_output(BM_output_file_list_with_path, combined_output):

    HGT_identity_dict = {}
    HGT_end_match_dict = {}
    HGT_full_length_match_dict = {}
    HGT_occurence_dict = {}
    HGT_concatenated_list = []

    for pwd_BM_output_file in BM_output_file_list_with_path:
        file_path, file_name = os.path.split(pwd_BM_output_file)
        taxon_rank = file_name.split('_')[-3][0]
        for BM_HGT in open(pwd_BM_output_file):
            if not BM_HGT.startswith('Gene_1'):
                BM_HGT_split = BM_HGT.strip().split('\t')

                gene_1 = BM_HGT_split[0]
                gene_2 = BM_HGT_split[1]
                identity = float(BM_HGT_split[4])
                end_match = BM_HGT_split[5]
                full_length_match = BM_HGT_split[6]
                concatenated = '%s___%s' % (gene_1, gene_2)

                if concatenated not in HGT_concatenated_list:
                    HGT_concatenated_list.append(concatenated)

                # store in dict
                if concatenated not in HGT_identity_dict:
                    HGT_identity_dict[concatenated] = identity

                if concatenated not in HGT_end_match_dict:
                    HGT_end_match_dict[concatenated] = end_match

                if concatenated not in HGT_full_length_match_dict:
                    HGT_full_length_match_dict[concatenated] = full_length_match


                if concatenated not in HGT_occurence_dict:
                    HGT_occurence_dict[concatenated] = [taxon_rank]
                else:
                    HGT_occurence_dict[concatenated].append(taxon_rank)


    detection_levels = ['p', 'c', 'o', 'f', 'g', 's']
    HGT_occurence_dict_formatted = {}
    for each_HGT in HGT_occurence_dict:
        occurence_str = ''
        for each_level in detection_levels:
            if each_level in HGT_occurence_dict[each_HGT]:
                occurence_str += '1'
            else:
                occurence_str += '0'
        HGT_occurence_dict_formatted[each_HGT] = occurence_str

    combined_output_handle = open(combined_output, 'w')
    combined_output_handle.write('Gene_1\tGene_2\tIdentity\toccurence(pcofgs)\tend_match\tfull_length_match\n')
    for concatenated_HGT in sorted(HGT_concatenated_list):
        concatenated_HGT_split = concatenated_HGT.split('___')


        for_out = '%s\t%s\t%s\t%s\t%s\t%s\n' % (concatenated_HGT_split[0],
                                                    concatenated_HGT_split[1],
                                                    HGT_identity_dict[concatenated_HGT],
                                                    HGT_occurence_dict_formatted[concatenated_HGT],
                                                    HGT_end_match_dict[concatenated_HGT],
                                                    HGT_full_length_match_dict[concatenated_HGT])
        combined_output_handle.write(for_out)

    combined_output_handle.close()


wd = '/Users/songweizhi/Desktop/555555/TT_90MGs'
BM_output_folder = 'TT_90MGs_BM'
combined_output = 'TT_90MGs_BM.txt'

output_prefix = 'GoodBins_0.5_0.05'
pwd_candidates_seq_file = 'GoodBins_0.5_0.05_all_combined_ffn.fasta'

pwd_candidates_file_ET_validated_fasta_nc = 'GoodBins_0.5_0.05_BM.ffn'
pwd_candidates_file_ET_validated_fasta_aa = 'GoodBins_0.5_0.05_BM.faa'


os.chdir(wd)


BM_output_file_re = '%s/*.txt' % BM_output_folder
BM_output_file_list_with_path = [('%s/%s' % (BM_output_folder, os.path.basename(file_name))) for file_name in glob.glob(BM_output_file_re)]


combine_BM_output(BM_output_file_list_with_path, combined_output)

validated_candidate_list = set()
for each in open(combined_output):
    if not each.startswith('Gene_1'):
        validated_candidate_list.add(each.strip().split('\t')[0])
        validated_candidate_list.add(each.strip().split('\t')[1])
print(len(validated_candidate_list))

# export sequence of validated candidates
combined_output_validated_fasta_nc_handle = open(pwd_candidates_file_ET_validated_fasta_nc, 'w')
combined_output_validated_fasta_aa_handle = open(pwd_candidates_file_ET_validated_fasta_aa, 'w')
for each_candidate in SeqIO.parse(pwd_candidates_seq_file, 'fasta'):
    if each_candidate.id in validated_candidate_list:
        # output nc sequences
        SeqIO.write(each_candidate, combined_output_validated_fasta_nc_handle, 'fasta')
        # output aa sequences
        each_candidate_aa = each_candidate
        each_candidate_aa.seq = each_candidate_aa.seq.translate()
        SeqIO.write(each_candidate_aa, combined_output_validated_fasta_aa_handle, 'fasta')
combined_output_validated_fasta_nc_handle.close()
combined_output_validated_fasta_aa_handle.close()


