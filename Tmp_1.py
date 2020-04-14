import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def get_bar_plot(num_list_1, num_list_2, label_list, label_rotation, output_plot):

    # Set position of bar on X axis
    barWidth = 0.4
    r1 = np.arange(len(num_list_1))
    r2 = [x + barWidth for x in r1]

    plt.bar(r1, num_list_1, color='orange', width=barWidth, edgecolor='white', label='Kelp')
    plt.bar(r2, num_list_2, color='lightblue', width=barWidth, edgecolor='white', label='Tara')
    plt.xticks([r + barWidth for r in range(len(num_list_1))], label_list, rotation=label_rotation, fontsize=5)
    plt.title('')
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.close()


# input files
cog_cate_to_compare =       'q'
cog_id_to_function =        '/Users/songweizhi/Desktop/diff_COG/combined.txt'
kelp_cog_file =             '/Users/songweizhi/Desktop/diff_COG/Kelp_%s.txt'                % cog_cate_to_compare
tara_cog_file =             '/Users/songweizhi/Desktop/diff_COG/Tara_%s.txt'                % cog_cate_to_compare

# output files
output_text_more_in_kelp =  '/Users/songweizhi/Desktop/diff_COG/COG_%s_more_in_kelp.txt'    % cog_cate_to_compare
output_text_more_in_tara =  '/Users/songweizhi/Desktop/diff_COG/COG_%s_more_in_tara.txt'    % cog_cate_to_compare
output_plot_more_in_kelp =  '/Users/songweizhi/Desktop/diff_COG/COG_%s_more_in_kelp.png'    % cog_cate_to_compare
output_plot_more_in_tara =  '/Users/songweizhi/Desktop/diff_COG/COG_%s_more_in_tara.png'    % cog_cate_to_compare


cog_id_to_function_dict = {}
for each_cog in open(cog_id_to_function):
    each_cog_split = each_cog.strip().split('\t')
    if (len(each_cog_split) > 1) and (each_cog_split[1] not in cog_id_to_function_dict):
        cog_id_to_function_dict[each_cog_split[1]] = each_cog_split[3]

uniq_both = []
cog_kelp_list = []
for cog_kelp in open(kelp_cog_file):
    cog_kelp_strip = cog_kelp.strip()
    cog_kelp_list.append(cog_kelp_strip)
    if cog_kelp_strip not in uniq_both:
        uniq_both.append(cog_kelp_strip)

cog_tara_list = []
for cog_tara in open(tara_cog_file):
    cog_tara_strip = cog_tara.strip()
    cog_tara_list.append(cog_tara_strip)
    if cog_tara_strip not in uniq_both:
        uniq_both.append(cog_tara_strip)

num_list_kelp = []
num_list_tara = []
uniq_both_sorted_more_in_kelp = []
num_list_kelp_more_in_kelp = []
num_list_tara_more_in_kelp = []
uniq_both_sorted_more_in_tara = []
num_list_kelp_more_in_tara = []
num_list_tara_more_in_tara = []
output_text_more_in_kelp_handle = open(output_text_more_in_kelp, 'w')
output_text_more_in_tara_handle = open(output_text_more_in_tara, 'w')
for uniq_cog in sorted(uniq_both):

    uniq_cog_num_pct_kelp = cog_kelp_list.count(uniq_cog)*100/len(cog_kelp_list)
    uniq_cog_num_pct_tara = cog_tara_list.count(uniq_cog)*100/len(cog_tara_list)
    uniq_cog_num_pct_kelp = float("{0:.3f}".format(uniq_cog_num_pct_kelp))
    uniq_cog_num_pct_tara = float("{0:.3f}".format(uniq_cog_num_pct_tara))
    num_list_kelp.append(uniq_cog_num_pct_kelp)
    num_list_tara.append(uniq_cog_num_pct_tara)

    if uniq_cog_num_pct_kelp > uniq_cog_num_pct_tara:
        if ((uniq_cog_num_pct_kelp - uniq_cog_num_pct_tara) / uniq_cog_num_pct_kelp >= 0.5) and (cog_kelp_list.count(uniq_cog) >=3):
            uniq_both_sorted_more_in_kelp.append(uniq_cog)
            num_list_kelp_more_in_kelp.append(uniq_cog_num_pct_kelp)
            num_list_tara_more_in_kelp.append(uniq_cog_num_pct_tara)
            output_text_more_in_kelp_handle.write('%s\t%s\t%s\t%s\n' % (uniq_cog, cog_id_to_function_dict[uniq_cog], uniq_cog_num_pct_kelp, uniq_cog_num_pct_tara))

    if uniq_cog_num_pct_tara > uniq_cog_num_pct_kelp:
        if ((uniq_cog_num_pct_tara - uniq_cog_num_pct_kelp) / uniq_cog_num_pct_tara >= 0.5) and (cog_tara_list.count(uniq_cog) >=3):
            uniq_both_sorted_more_in_tara.append(uniq_cog)
            num_list_kelp_more_in_tara.append(uniq_cog_num_pct_kelp)
            num_list_tara_more_in_tara.append(uniq_cog_num_pct_tara)
            output_text_more_in_tara_handle.write('%s\t%s\t%s\t%s\n' % (uniq_cog, cog_id_to_function_dict[uniq_cog], uniq_cog_num_pct_kelp, uniq_cog_num_pct_tara))

output_text_more_in_kelp_handle.close()
output_text_more_in_tara_handle.close()

get_bar_plot(num_list_kelp_more_in_kelp, num_list_tara_more_in_kelp, uniq_both_sorted_more_in_kelp, 270, output_plot_more_in_kelp)
get_bar_plot(num_list_kelp_more_in_tara, num_list_tara_more_in_tara, uniq_both_sorted_more_in_tara, 270, output_plot_more_in_tara)
