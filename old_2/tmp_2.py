


done_list = []
for each in open('/Users/songweizhi/Desktop/000.txt'):
    done_list.append(each.strip())


all_list = []
for each in open('/Users/songweizhi/Desktop/001.txt'):
    all_list.append(each.strip())

for each in all_list:

    if each not in done_list:

        print('cp Tara_SD_all_prodigal_output/%s.faa Tara_SD_all_prodigal_output_rest/ &' % each)