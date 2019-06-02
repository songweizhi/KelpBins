

for each in open('/Users/songweizhi/Desktop/mapping_cmds.sh'):

    sam_file_id = each.strip().split(' -S ')[1].split('.sam -p')[0]
    file_name = '/Users/songweizhi/Desktop/qsub_mapping/qsub_mapping_%s.sh' % sam_file_id

    file_name_handle = open(file_name, 'w')
    file_name_handle.write('#!/bin/bash\n#PBS -l nodes=1:ppn=6\n#PBS -l mem=60gb\n#PBS -l walltime=23:59:00\n#PBS -j oe\n#PBS -M songwz03@163.com\n#PBS -m ae\n\nmodule load bowtie/2.3.4.2\nmodule load samtools/1.9\n\n')
    file_name_handle.write('cd /srv/scratch/z5039045/Liu_RL/Mariana_Trench/04_mapping\n')
    file_name_handle.write(each)
    file_name_handle.write('samtools view -bS %s.sam -o %s.bam\n' % (sam_file_id, sam_file_id))
    file_name_handle.write('rm %s.sam\n' % sam_file_id)
    file_name_handle.write('samtools sort %s.bam -o %s_sorted.bam\n' % (sam_file_id, sam_file_id))
    file_name_handle.write('rm %s.bam\n' % sam_file_id)
    file_name_handle.write('samtools index %s_sorted.bam\n' % sam_file_id)
    file_name_handle.close()
