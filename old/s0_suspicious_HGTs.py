import os
import glob
from s0_Kelp_bins_config import suspicious_HGTs_plot_folder


suspicious_HGTs_file_re = '%s/*.eps' % suspicious_HGTs_plot_folder
suspicious_HGTs_file_list = [os.path.basename(file_name) for file_name in glob.glob(suspicious_HGTs_file_re)]

suspicious_HGTs = ['.'.join(i.split('.')[:-1]) for i in suspicious_HGTs_file_list]

#print(suspicious_HGTs)