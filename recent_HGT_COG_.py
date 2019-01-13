
HGT_PG_validated_o =         '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05_o34_HGTs_PG_validated.txt'
PG_validated_o_gene_id_cog = '/Users/songweizhi/Desktop/KelpBins/COG_enrichment_plot_2/GoodBins_0.5_0.05_o34_HGTs_PG_aa_COG_results/protein-id_cog.txt'
PG_validated_o_cog_stats =   '/Users/songweizhi/Desktop/KelpBins/COG_enrichment_plot_2/GoodBins_0.5_0.05_o34_HGTs_PG_aa_COG_results/cog_stats.txt'


gene_to_cog_dict = {}
for each_gene_COG in open(PG_validated_o_gene_id_cog):
    each_gene_COG_split = each_gene_COG.strip().split('\t')
    gene_id= each_gene_COG_split[0]
    COG_id= each_gene_COG_split[1]
    gene_to_cog_dict[gene_id] = COG_id


cog_id_to_description_dict = {}
for each_cog_id in open(PG_validated_o_cog_stats):
    each_cog_id_split = each_cog_id.strip().split('\t')
    cog_id = each_cog_id_split[0]
    cog_description = each_cog_id_split[1]
    cog_id_to_description_dict[cog_id] = cog_description


for hgt_o in open(HGT_PG_validated_o):
    if not hgt_o.startswith('Gene_1'):
        hgt_o_split = hgt_o.strip().split('\t')
        gene_1 = hgt_o_split[0]
        gene_2 = hgt_o_split[1]
        identity = float(hgt_o_split[4])

        gene_1_cog = 'NA'
        gene_2_cog = 'NA'
        gene_1_description = 'NA'
        gene_2_description = 'NA'

        if gene_1 in gene_to_cog_dict:
            gene_1_cog = gene_to_cog_dict[gene_1]
        if gene_2 in gene_to_cog_dict:
            gene_2_cog = gene_to_cog_dict[gene_2]
        if gene_1_cog in cog_id_to_description_dict:
            gene_1_description = cog_id_to_description_dict[gene_1_cog]
        if gene_2_cog in cog_id_to_description_dict:
            gene_2_description = cog_id_to_description_dict[gene_2_cog]

        for_out_gene_1 = '%s\t%s\t%s' % (gene_1, gene_1_cog, gene_1_description)
        for_out_gene_2 = '%s\t%s\t%s' % (gene_2, gene_2_cog, gene_2_description)

        if identity >= 90:
            print(for_out_gene_1)
            print(for_out_gene_2)
