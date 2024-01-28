def search_seq(output_file, homolog_file, chrom, type_del ):
    # Define a mapping of chromosome names to numerical values
    chr_list = ['1','2','3','4', '5', '6', '7', '8', '9','10', '11', '12', 
                '13', '14', '15', '16', '17', '18','19','20', '21',  '22', 'X', 'Y']
    num_list = list(range(0,24))
    chr_key = dict(zip(chr_list, num_list))

    # Path to the reference genome in FASTA format
    ref_genome_chr = '/home/labs/shlush/shared/broad_hg19/Homo_sapiens_assembly19.fa'
    
   
    chr_seq = list(SeqIO.parse(ref_genome_chr, "fasta"))
    
    # List to store generated homolog sequences
    count_total = []  
    
   
    ofile = open(homolog_file, 'w')
    
    # Read homolog information from the input file
    with open(homolog_file, 'r') as f:
        f.readline()  # Skip the header line
        for i, line in enumerate(f):
            # Parse homolog information (Start, End, Homolog)
            Start_del, End_del, Homolog = (int(x) for x in line.strip().split())
    
           
            Chr = chr_key[chrom]
            
            # Define genomic regions for homolog sequence extraction

            Start_seq1 = Start_del - 15
            End_seq1 = End_del + 15
            
            # Extract left and right sequences based on deletion type
            if type_del == 'canonical':
                str1_left = chr_seq[Chr][int(Start_seq1-1):int(Start_del)].seq 
                str1_right = chr_seq[Chr][int(End_del):int(End_seq1)].seq
            elif type_del == 'left':
                str1_left = chr_seq[Chr][int(Start_seq1-1):int(Start_del-1)].seq 
                str1_right = chr_seq[Chr][int(End_del):int(End_seq1)].seq
            elif type_del == 'right':
                str1_left = chr_seq[Chr][int(Start_seq1-1):int(Start_del)+int(Homolog)].seq 
                str1_right = chr_seq[Chr][int(End_del)+int(Homolog+1):int(End_seq1)].seq
    
            # Concatenate left and right flanks of the deletion
            str1 = str(str1_left) + str(str1_right)
            
            
            count_total.append('>{CHR}_{s}_{e}_{l}\n'.format(CHR=chrom, s=Start_del, e=End_del, l=Homolog) + str1 + "\n")
    
            if i % 100000 == 0:
                ofile.writelines(count_total)
                count_total = []
    
        # Write any remaining sequences to the output file
        ofile.writelines(count_total)
        ofile.close()

