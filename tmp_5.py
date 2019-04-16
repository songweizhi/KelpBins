import s0_Kelp_bins_config


direction_1_way = 0
direction_2_ways = 0


for each_HGT in open(s0_Kelp_bins_config.HGT_PG_validated_txt_pcofg):

    if not each_HGT.startswith('Gene_1'):
        each_HGT_split = each_HGT.strip().split('\t')

        occurence = each_HGT_split[3]
        direction = each_HGT_split[6]


        # print(occurence)
        # print(direction)
        if occurence.count('1') > 1:
            print('%s\t%s\t%s\t%s' % (each_HGT_split[0], each_HGT_split[1], occurence, direction))

            if direction == 'both':
                direction_2_ways += 1
            else:
                direction_1_way += 1


    # print(each_HGT)
print(direction_1_way)
print(direction_2_ways)
