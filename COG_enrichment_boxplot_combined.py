import os
import glob
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


# file in
COG_func_stats_folder =     '/Users/songweizhi/Desktop/Kelp_NM/COG/COG_func_stats_combined'
Kelp_HGT_func_stats_file =  '/Users/songweizhi/Desktop/Kelp_NM/COG/zKelp_recipient_genes_func_stats_GeneNumber.txt'
Tara_HGT_func_stats_file =  '/Users/songweizhi/Desktop/Kelp_NM/COG/zTara_NM_recipient_genes_func_stats_GeneNumber.txt'
column_order = ['J', 'K', 'L', 'V', 'T', 'M', 'U', 'O', 'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S']

# file out
df_out = '/Users/songweizhi/Desktop/Kelp_NM/COG/COG_enrichment_boxplot_combined.txt'
png_out = '/Users/songweizhi/Desktop/Kelp_NM/COG/COG_enrichment_boxplot_combined.png'


COG_func_stats_file_re = '%s/*.txt' % COG_func_stats_folder
COG_func_stats_file_list = [os.path.basename(file_name) for file_name in glob.glob(COG_func_stats_file_re)]

df_out_handle = open(df_out, 'w')
df_out_handle.write('MAG,COG_cate,Percent,Lifestyle\n')
for Kelp_COG_func_stats in COG_func_stats_file_list:

    pwd_Kelp_COG_func_stats = '%s/%s' % (COG_func_stats_folder, Kelp_COG_func_stats)
    MAG_id = '_'.join(Kelp_COG_func_stats.split('_')[:-3])

    MAG_source = 'Planktonic'
    if 'Refined' in MAG_id:
        MAG_source = 'Kelp-associated'

    current_genome_cate_num_dict = {}
    current_genome_cate_total_num = 0
    for cate in open(pwd_Kelp_COG_func_stats):
        if not cate.startswith('Category'):
            cate_split = cate.strip().split('\t')
            current_genome_cate_num_dict[cate_split[0]] = int(cate_split[1])
            current_genome_cate_total_num += int(cate_split[1])

    current_genome_cate_pct_dict = {}
    for i in current_genome_cate_num_dict:
        i_pct = float("{0:.2f}".format(current_genome_cate_num_dict[i]*100/current_genome_cate_total_num))
        current_genome_cate_pct_dict[i] = i_pct

    for cate_pct in current_genome_cate_pct_dict:
        if cate_pct in column_order:
            df_out_handle.write('%s,%s,%s,%s\n' % (MAG_id, cate_pct, str(current_genome_cate_pct_dict[cate_pct]), MAG_source))
df_out_handle.close()


# get percentage among HGT
Kelp_HGT_cate_num_dict = {}
Kelp_HGT_cate_total_num = 0
for Kelp_cate in open(Kelp_HGT_func_stats_file):
    if not Kelp_cate.startswith('Category'):
        Kelp_cate_split = Kelp_cate.strip().split('\t')
        Kelp_HGT_cate_num_dict[Kelp_cate_split[0]] = int(Kelp_cate_split[1])
        Kelp_HGT_cate_total_num += int(Kelp_cate_split[1])

Tara_HGT_cate_num_dict = {}
Tara_HGT_cate_total_num = 0
for Tara_cate in open(Tara_HGT_func_stats_file):
    if not Tara_cate.startswith('Category'):
        Tara_cate_split = Tara_cate.strip().split('\t')
        Tara_HGT_cate_num_dict[Tara_cate_split[0]] = int(Tara_cate_split[1])
        Tara_HGT_cate_total_num += int(Tara_cate_split[1])

Kelp_HGT_cate_pct_dict = {}
for i in Kelp_HGT_cate_num_dict:
    i_pct = float("{0:.2f}".format(Kelp_HGT_cate_num_dict[i]*100/Kelp_HGT_cate_total_num))
    Kelp_HGT_cate_pct_dict[i] = i_pct

Tara_HGT_cate_pct_dict = {}
for j in Tara_HGT_cate_num_dict:
    j_pct = float("{0:.2f}".format(Tara_HGT_cate_num_dict[j]*100/Tara_HGT_cate_total_num))
    Tara_HGT_cate_pct_dict[j] = j_pct


Kelp_HGT_cate_pct_color_list = ['red', 'deepskyblue', 'deepskyblue', 'grey', 'deepskyblue', 'deepskyblue', 'deepskyblue', 'deepskyblue', 'red', 'red', 'red', 'grey', 'grey', 'red', 'red', 'red', 'grey', 'deepskyblue']
Tara_HGT_cate_pct_color_list = ['red', 'grey', 'deepskyblue', 'grey', 'grey', 'deepskyblue', 'deepskyblue', 'deepskyblue', 'red', 'grey', 'red', 'grey', 'grey', 'grey', 'red', 'grey', 'grey', 'deepskyblue']


Kelp_HGT_cate_pct_shape_list = []
for color in Kelp_HGT_cate_pct_color_list:
    if color == 'red':
        Kelp_HGT_cate_pct_shape_list.append('^')
    if color == 'black':
        Kelp_HGT_cate_pct_shape_list.append('s')
    if color == 'deepskyblue':
        Kelp_HGT_cate_pct_shape_list.append('v')
    if color == 'grey':
        Kelp_HGT_cate_pct_shape_list.append('s')


Tara_HGT_cate_pct_shape_list = []
for color in Tara_HGT_cate_pct_color_list:
    if color == 'red':
        Tara_HGT_cate_pct_shape_list.append('^')
    if color == 'black':
        Tara_HGT_cate_pct_shape_list.append('s')
    if color == 'deepskyblue':
        Tara_HGT_cate_pct_shape_list.append('v')
    if color == 'grey':
        Tara_HGT_cate_pct_shape_list.append('s')


MAG_df = pd.read_csv(df_out)

box_order =   ['Kelp-associated',   'Planktonic']
color_order = ['orange',            'lightblue']

sns.boxplot(data=MAG_df, x="COG_cate", y="Percent", order=column_order,
            hue="Lifestyle", hue_order=box_order, palette=color_order, fliersize=1.2, showfliers=0)

# add dots for kelp HGTs
n = 0
for cate in column_order:
    plt.plot(n - 0.2, Kelp_HGT_cate_pct_dict[cate], alpha=1, marker=Kelp_HGT_cate_pct_shape_list[n], markersize=5, markeredgewidth=0, color=Kelp_HGT_cate_pct_color_list[n])
    n += 1

# add dots for Tara HGTs
m = 0
for cate in column_order:
    plt.plot(m + 0.2, Tara_HGT_cate_pct_dict[cate], alpha=1, marker=Tara_HGT_cate_pct_shape_list[m], markersize=5, markeredgewidth=0, color=Tara_HGT_cate_pct_color_list[m])
    m += 1

plt.tight_layout()
plt.savefig(png_out, bbox_inches='tight', dpi=300)
plt.close()
