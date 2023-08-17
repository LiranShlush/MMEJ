#This code takes the MMEJ search sequences as an input and outputs the number of mutant reads and the coverage in that region for each MMEJ deletion event 


# The input_file_name is a file containing details related to the cooridnates of the homolog
# The following columns are required:

# Column 1 – Chromsome(as given in the reference genome)
# Column 2 - Start of the homolog1
# Column 3 – Start of the homolog2
# Column 4 – Length of the homolog
# Column 5 - Canonical MMEJ search sequence
# Column 6 – Imperfect type A MMEJ search sequence
# Column 7 – Imperfect type B MMEJ search sequence


input_file = str(sys.argv[1])
bam_name = str(sys.argv[2])
output_file = str(sys.argv[3])


out = open(output_file, 'w')
count_total = []
with pysam.AlignmentFile(bam_name, "rb") as outf:
    with open(input_file,'r') as f:
        p = f.readline()
        var_names = p.strip().split('\t')

        for i,line in enumerate(f):
            row = dict(zip(var_names,  line.strip().split('\t')))

            
            Start_del = int(float(row['Start1']))- 1
            End_del = int(float(row['Start2']))

            Start_del1 = Start_del-100
            End_del1 = End_del + 100
            str_canonical =  str(row.get('canonical'))
            str_left =  str(row.get('left'))
            str_right =  str(row.get('right'))
            #For hg18
            # Chr = "chr" + str(row.get('Chr'))
            #for hg19
            Chr =  str(row.get('Chr'))

            start_5 = outf.count_coverage(Chr, start=Start_del-5, end=Start_del-4)
            end_5 = outf.count_coverage(Chr.format(Chr = Chr), start=End_del+5, end=End_del+6)
            basecov_start_5 = start_5[0][0] + start_5[1][0] + start_5[2][0] + start_5[3][0]
            basecov_end_5 = end_5[0][0] + end_5[1][0] + end_5[2][0] + end_5[3][0]
            mean_coverage = round(np.mean([basecov_start_5, basecov_end_5]))
            
            count_canonical = 0
            count_left = 0
            count_right = 0
            
            try:

                x =  outf.fetch(Chr, Start_del1, End_del1)
                for read in x:

                    if (str_canonical  != ""):
                        
                        if (str_canonical in read.query_sequence):
                            count_canonical = count_canonical + 1

                    if (str_left  != ""):

                        if (str_left in read.query_sequence):
                            count_left = count_left + 1

                    if (count_right  != ""):
                        if (str_right in read.query_sequence):
                            count_right = count_right + 1
     
                   
                    
            except KeyboardInterrupt:
                print("An TypeError occurred:", i) 

            else:
                pass

            count_total.append('{i}\t{c}\t{S}\t{E}\t{L}\t{count_canonical}\t{count_left}\t{count_right}\t{mean}\n'.format(i=i,
                                          c = row.get('Chr'),   S = row['Start'], E = row['End'], L = row['Length'], \
                                             mean = mean_coverage,  count_canonical=count_canonical, \
                                         count_left =count_left,  count_right = count_right ))
            if i % 100000:


                out.writelines(count_total)
                count_total = []

        outf.close()
        out.writelines(count_total)
        out.close()

