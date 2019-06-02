import os
import glob
from Bio import SeqIO
from copy import deepcopy


bin_folder =                    '/Users/songweizhi/Desktop/KelpBins/GoodBins_0.5_0.05'


bin_file_re = '%s/*.fasta' % bin_folder

bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]



current_clustet_combined_genome = '/Users/songweizhi/Desktop/combined_528bins.fasta'


current_clustet_combined_genome_handle = open(current_clustet_combined_genome, 'w')
for each_genome in bin_file_list:
    genome_name_no_ext = '.'.join(each_genome.split('.')[:-1])

    # combine current cluster genomes and add genome name to contig ID
    pwd_genome = '%s/%s' % (bin_folder, each_genome)
    for contig in SeqIO.parse(pwd_genome, 'fasta'):
        contig_new = deepcopy(contig)
        contig_new.id = '%s___%s' % (genome_name_no_ext, str(contig.id))
        SeqIO.write(contig_new, current_clustet_combined_genome_handle, 'fasta')

current_clustet_combined_genome_handle.close()

