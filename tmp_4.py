import scipy
from scipy import stats
import pandas as pd
from pandas import DataFrame
from matplotlib import pyplot


group1 = [1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,8,8,9,0, 10,10,10,10,10, 10,10,10,10,10, 10, 10, 10]

group1_normality = 'No'
if float(stats.normaltest(group1).pvalue) > 0.05:
    group1_normality = 'Yes'


print('Normality Test: %s' % group1_normality)
pyplot.hist(group1)
pyplot.show()