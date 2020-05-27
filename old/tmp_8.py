
protein_id_cog_file = '/Users/songweizhi/Desktop/ResistFlow_pcofg_detected_HGTs_recipient_genes_COG_annot/ResistFlow_pcofg_detected_HGTs_recipient_genes_protein-id_cog.txt'
cog_stats_file = '/Users/songweizhi/Desktop/ResistFlow_pcofg_detected_HGTs_recipient_genes_COG_annot/ResistFlow_pcofg_detected_HGTs_recipient_genes_cog_stats.txt'


cog_id_to_description_dict = {}
for cog_id in open(cog_stats_file):
    cog_id_split = cog_id.strip().split('\t')
    cog_id_to_description_dict[cog_id_split[0]] = cog_id_split[1]


out_handle = open('/Users/songweizhi/Desktop/gene_to_cog.txt', 'w')
cog_cate_to_id_dict = {}
for protein_to_cog in open(protein_id_cog_file):
    protein_id, cog_id, *cog_cates = protein_to_cog.strip().split('\t')
    out_handle.write('%s\t%s\t%s\t%s\n' % (protein_id, cog_id, ''.join(cog_cates), cog_id_to_description_dict[cog_id]))
out_handle.close()
