import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


################################################# 10 Minutes to pandas #################################################

# https://pandas.pydata.org/pandas-docs/stable/10min.html#min

# Creating a Series by passing a list of values, letting pandas create a default integer index
s = pd.Series([1, 3, 5, np.nan, 6, 8])
print(s)
print()


# Creating a DataFrame by passing a NumPy array, with a datetime index and labeled columns
dates = pd.date_range('20130101', periods=6, freq='D')
print(dates)
print()


# Creating a DataFrame by passing a NumPy array, with a datetime index and labeled columns
num_list = [[0, 25, 33, 42],
            [41, 22, 23, 46],
            [32,42, 23, 49],
            [21, 12, 43, 40],
            [15, 12, 31, 34],
            [3, 2, 32, 64]]

df = pd.DataFrame(np.array(num_list), index=dates, columns=list('ABCD'))
print(df)
print()

print(df)

# Creating a DataFrame by passing a dict of objects that can be converted to series-like
df2 = pd.DataFrame({'A': 1.,
                    'B': pd.Timestamp('20130102'),
                    'C': pd.Series(1, index=list(range(4)), dtype='float32'),
                    'D': np.array([3] * 4, dtype='int32'),
                    'E': pd.Categorical(["test", "train", "test", "train"]),
                    'F': 'foo'})

print(df2)
print()
print(df2.dtypes)
print()


# view the top and bottom rows of the DataFrame
print(df.head(2))
print(df.tail(3))
print()

# Display the index, columns, and the underlying NumPy data
print(df.index)
print()
print(df.columns)
print()
print(df.values)
print()


# show a quick statistic summary of your data:
print(df.describe())
print()


# Transposing your data
df_t = df.T
print(df_t)
print()


# Sorting by an axis (column name)
print(df.sort_index(axis=1, ascending=False))
print()


# Sorting by values
print(df.sort_values(by='B'))
print()


# Selecting a single column, equivalent to df.A
print(df['A'])
print()
print(df.A)
print()


# Selecting via [], which slices the rows
print(df[0:3])
print()
print(df['20130102':'20130104'])
print()


# Selection by Label
print(df.loc[dates[0]])
print()


# Selecting on a multi-axis by label
print(df.loc[:, ['A', 'C']])
print()


# Showing label slicing, both endpoints are included
print(df.loc['20130102':'20130104', ['A', 'C']])
print()


# Reduction in the dimensions of the returned object
print(df.loc['20130102',['A', 'B']])
print()


# For getting a scalar value
print(df.loc[dates[0], 'A'])
print()


# For getting fast access to a scalar (equivalent to the prior method)
print(df.at[dates[0],'A'])
print()


print(df)
print()


df_D = df.loc[:, df.columns == 'D']
print(df_D)
print()

print(df_D.iloc[:,0])
print(df_D['D'].tolist())



print()


print(df)

df2 = df.drop(['B', 'C'], axis=1)

print(df2)





#print(df.loc[:, df.columns == 'D'].values)



# print(df.loc[:, df.columns != 'D'])

# df.boxplot(vert=False, grid=False, figsize=(20, 10))
# plt.savefig('/Users/songweizhi/Desktop/test2.png', dpi=300)
#

#################################################### Pandas Cookbook ###################################################

# https://pandas.pydata.org/pandas-docs/stable/cookbook.html#cookbook


# df = pd.DataFrame({'AAA': [4, 5, 6, 7],
#                    'BBB': [10, 20, 30, 40],
#                    'CCC': [100, 50, -30, -50]})
# print(df)
# print()
#
#
# df.loc[df.AAA >= 5, 'BBB'] = -1
# print(df)
# print()
#
#
# df.loc[df.AAA >= 5,['BBB', 'CCC']] = 555
# print(df)
# print()
#
#
# df.loc[df.AAA < 5,['BBB', 'CCC']] = 2000
# print(df)
# print()
#
#
# df_mask = pd.DataFrame({'AAA': [True] * 4, 'BBB': [False] * 4, 'CCC': [True, False] * 2})
# print(df_mask)
# print()
# print(df.where(df_mask, -1000))
# print()

