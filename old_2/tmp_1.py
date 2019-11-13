

Kelp_all_BM_candidate_list = []
for each in open('/Users/songweizhi/Desktop/333333/Kelp_c15_HGTs_BM.txt'):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        Kelp_all_BM_candidate_list.append(each_split[0])


Kelp_dRep_BM_candidate_list = []
for each in open('/Users/songweizhi/Desktop/333333/Kelp_dRep_c15_HGTs_BM.txt'):
    if not each.startswith('Gene_1'):
        each_split = each.strip().split('\t')
        Kelp_dRep_BM_candidate_list.append(each_split[0])


shared_HGTs = set(Kelp_all_BM_candidate_list).intersection(Kelp_dRep_BM_candidate_list)


Kelp_all_uniq_HGT_list = []
for HGT in Kelp_all_BM_candidate_list:    
    if HGT not in shared_HGTs:
        Kelp_all_uniq_HGT_list.append(HGT)
        
        
Kelp_dRep_uniq_HGT_list = []
for HGT in Kelp_dRep_BM_candidate_list:
    if HGT not in shared_HGTs:
        Kelp_dRep_uniq_HGT_list.append(HGT)


print('Kelp_all_uniq_HGT_list')
print(Kelp_all_uniq_HGT_list)
print()
print('Kelp_dRep_uniq_HGT_list')
print(Kelp_dRep_uniq_HGT_list)
print()

target_gene = 'CB_ER_130617_Refined_4_02807'


for gene in open('/Users/songweizhi/Desktop/333333/Kelp_c15_subjects_in_one_line.txt'):
    gene_split = gene.strip().split('\t')
    gene_id = gene_split[0].split('|')[1]
    if gene_id == target_gene:
        subject_renamed_list = []
        for subject in gene_split[1:]:
            subject_renamed = '%s|%s' % (subject.split('|')[0], subject.split('|')[2])
            subject_renamed_list.append(subject_renamed)
        print(gene.strip())
        print('All\t%s\t%s' % (gene_split[0], '\t'.join(subject_renamed_list)))

for gene in open('/Users/songweizhi/Desktop/333333/Kelp_dRep_c15_subjects_in_one_line.txt'):
    gene_split = gene.strip().split('\t')
    gene_id = gene_split[0].split('|')[1]
    if gene_id == target_gene:
        subject_renamed_list = []
        for subject in gene_split[1:]:
            subject_renamed = '%s|%s' % (subject.split('|')[0], subject.split('|')[2])
            subject_renamed_list.append(subject_renamed)
        print('dRep\t%s\t%s' % (gene_split[0], '\t'.join(subject_renamed_list)))


# CB_ER_130617_Refined_4_02807
# BH_ER_161216_Refined_12_01674

print(Kelp_all_BM_candidate_list.count('CB_ER_130617_Refined_4_02807'))
print(Kelp_all_BM_candidate_list.count('BH_ER_161216_Refined_12_01674'))


