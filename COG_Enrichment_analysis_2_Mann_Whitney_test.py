import scipy
from scipy import stats
import pandas as pd
from pandas import DataFrame
from matplotlib import pyplot
from statsmodels.stats.multitest import multipletests


def remove_0_from_Pandas_Series(Pandas_Series):

    no_0_num_list = []
    for index, value in Pandas_Series.items():
        if value > 0:
            no_0_num_list.append(value)

    return no_0_num_list


# file in
#csv_file =     '/Users/songweizhi/Desktop/COG_enrichment_analysis/Herbivore_specialisation.csv'
csv_file =      '/Users/songweizhi/Desktop/COG_enrichment_analysis/COG_enrichment_analysis.csv'
output_test =   '/Users/songweizhi/Desktop/COG_enrichment_analysis/COG_enrichment_test_results_Mann_Whitney_U.tab'


# df = pd.read_csv(csv_file, index_col=0)
df = pd.read_csv(csv_file)

# get a list of all columns in the dataframe without the Group column
column_list = [x for x in df.columns if x != 'Group']

# loop over column_list and execute code explained above
Mann_Whitney_U_test_results = {}
cog_to_group_mean_dict = {}
cog_to_group_mean_dict_no_zero = {}
n = 0
p_value_list = []
COG_id_list = []
for column in column_list:

    if n >= 2:

        group1 = df.where(df.Source == 'kelp-associated').dropna()[column]
        group2 = df.where(df.Source == 'planktonic').dropna()[column]
        #group1_no_zero = remove_0_from_Pandas_Series(group1)
        #group2_no_zero = remove_0_from_Pandas_Series(group2)

        #print('Processing the %sth column: %s' % (n, column))

        # store group mean into dict
        group1_mean = float("{0:.2f}".format(sum(group1) / len(group1)))
        group2_mean = float("{0:.2f}".format(sum(group2) / len(group2)))
        #group1_no_zero_mean = float("{0:.2f}".format(sum(group1_no_zero) / len(group1_no_zero)))
        #group2_no_zero_mean = float("{0:.2f}".format(sum(group2_no_zero) / len(group2_no_zero)))
        cog_to_group_mean_dict[column] = [group1_mean, group2_mean]
        #cog_to_group_mean_dict_no_zero[column] = [group1_no_zero_mean, group2_no_zero_mean]

        # perform Mann-Whitney U test
        Mann_Whitney_U_test_results[column] = stats.mannwhitneyu(group1, group2)[1]

        COG_id_list.append(column)
        p_value_list.append(stats.ttest_ind(group1, group2, equal_var=False)[1])

    n += 1

p_value_list_adjusted = multipletests(p_value_list, alpha=0.1, method='bonferroni')[1]


x = 0
output_test_handle = open(output_test, 'w')
output_test_handle.write('COG\tKelp\tPlanktonic\tP_value\tP_value_adjusted\n')
while x < len(COG_id_list):

    current_p            = float("{0:.3f}".format(p_value_list[x]))
    current_p_adjusted   = float("{0:.3f}".format(p_value_list_adjusted[x]))
    current_p_group_mean = [str(i) for i in cog_to_group_mean_dict[COG_id_list[x]]]

    output_test_handle.write('%s\t%s\t%s\t%s\n' % (COG_id_list[x], '\t'.join(current_p_group_mean), current_p, current_p_adjusted))

    x += 1

output_test_handle.close()

