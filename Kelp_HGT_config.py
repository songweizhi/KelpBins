from Bio import SeqIO


########################################################################################################################

# wd
Kelp_NM_dRep99_wd =                         '/Users/songweizhi/Desktop/Kelp_NM'

# files in
Kelp_NM_dRep99_pcofg_detected_HGTs =     '%s/0_file_in/Kelp_NM_dRep99_pcofg_detected_HGTs.txt'                 % Kelp_NM_dRep99_wd
Kelp_NM_dRep99_pcofg_detected_HGT_seqs = '%s/0_file_in/Kelp_NM_dRep99_pcofg_detected_HGTs_recipient_genes.faa' % Kelp_NM_dRep99_wd
Kelp_NM_dRep99_genome_size_file =        '%s/0_file_in/Kelp_NM_dRep99_all_genome_size.txt'                     % Kelp_NM_dRep99_wd
MAG_source_file =                        '%s/0_file_in/MAG_source.txt'                                         % Kelp_NM_dRep99_wd

# files out
Kelp_recipient_gene_seq_file =           '%s/zKelp_recipient_genes.faa'                                         % Kelp_NM_dRep99_wd
Tara_NM_recipient_gene_seq_file =        '%s/zTara_NM_recipient_genes.faa'                                      % Kelp_NM_dRep99_wd


# script
Boxplot_last1row = '~/PycharmProjects/MyBioTools/MyBioTools/COG_boxplot_last1row.R'


########################################################################################################################

# functional annotation files in
Kelp_dbCAN_df_txt =     '%s/0_file_in/faa_files_Kelp_and_HGT_dbCAN_percentage.txt'                % Kelp_NM_dRep99_wd
Kelp_dbCAN_df_txt_t =   '%s/0_file_in/faa_files_Kelp_and_HGT_dbCAN_percentage_t.txt'              % Kelp_NM_dRep99_wd
Tara_dbCAN_df_txt =     '%s/0_file_in/faa_files_Tara_NM_and_HGT_dbCAN_percentage.txt'             % Kelp_NM_dRep99_wd
Tara_dbCAN_df_txt_t =   '%s/0_file_in/faa_files_Tara_NM_and_HGT_dbCAN_percentage_t.txt'           % Kelp_NM_dRep99_wd

Kelp_COG_df_txt =       '%s/0_file_in/faa_files_Kelp_and_HGT_cog_cate_percentage.txt'             % Kelp_NM_dRep99_wd
Kelp_COG_df_txt_t =     '%s/0_file_in/faa_files_Kelp_and_HGT_cog_cate_percentage_t.txt'           % Kelp_NM_dRep99_wd
Tara_COG_df_txt =       '%s/0_file_in/faa_files_Tara_NM_and_HGT_cog_cate_percentage.txt'          % Kelp_NM_dRep99_wd
Tara_COG_df_txt_t =     '%s/0_file_in/faa_files_Tara_NM_and_HGT_cog_cate_percentage_t.txt'        % Kelp_NM_dRep99_wd

Kelp_KEGG_df_txt =      '%s/0_file_in/faa_files_Kelp_and_HGT_KEGG_C_percentage.txt'               % Kelp_NM_dRep99_wd
Kelp_KEGG_df_txt_t =    '%s/0_file_in/faa_files_Kelp_and_HGT_KEGG_C_percentage_t.txt'             % Kelp_NM_dRep99_wd
Tara_KEGG_df_txt =      '%s/0_file_in/faa_files_Tara_NM_and_HGT_KEGG_C_percentage.txt'            % Kelp_NM_dRep99_wd
Tara_KEGG_df_txt_t =    '%s/0_file_in/faa_files_Tara_NM_and_HGT_KEGG_C_percentage_t.txt'          % Kelp_NM_dRep99_wd


# functional annotation files out
Kelp_dbCAN_df_png = '%s/faa_files_Kelp_and_HGT_dbCAN_percentage.png'        % Kelp_NM_dRep99_wd
Tara_dbCAN_df_png = '%s/faa_files_Tara_NM_and_HGT_dbCAN_percentage.png'     % Kelp_NM_dRep99_wd


########################################################################################################################

# get Kelp and Tara MAG list
Kelp_bin_list = []
Tara_NM_bin_list = []
for MAG in open(MAG_source_file):
    MAG_split = MAG.strip().split('\t')
    if MAG_split[0] == 'Kelp':
        Kelp_bin_list.append(MAG_split[1])
    if MAG_split[0] == 'Tara_NM':
        Tara_NM_bin_list.append(MAG_split[1])


# get the sequences of recipient genes for Kelp and Tara_NM MAGs
Kelp_recipient_gene_seq_file_handle = open(Kelp_recipient_gene_seq_file, 'w')
Tara_NM_recipient_gene_seq_file_handle = open(Tara_NM_recipient_gene_seq_file, 'w')
for recipient in SeqIO.parse(Kelp_NM_dRep99_pcofg_detected_HGT_seqs, 'fasta'):
    recipient_MAG = '_'.join(recipient.id.split('_')[:-1])
    if recipient_MAG in Kelp_bin_list:
        SeqIO.write(recipient, Kelp_recipient_gene_seq_file_handle, 'fasta')
    if recipient_MAG in Tara_NM_bin_list:
        SeqIO.write(recipient, Tara_NM_recipient_gene_seq_file_handle, 'fasta')
Kelp_recipient_gene_seq_file_handle.close()
Tara_NM_recipient_gene_seq_file_handle.close()
