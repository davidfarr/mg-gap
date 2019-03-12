# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 10:43:53 2019

@author: Heathro

Description: Reduces a vcf file to meta section and 
one line for each chromosome number for testing and
debugging purposes.
"""

# Open files to read from and write to
vcfpath = open("D:/MG_GAP/Ali_w_767.vcf", "rU")
testvcf = open("REDUCED_ali.vcf", "w")

# Keep track of chromosome number so we can get one of each
temp_chrom = 0
counter = 0

for line_index, line in enumerate(vcfpath):

    # Found a chromosome line
    if line[0:8] == "sNNffold":

        column = line.split('\t')
        first_col = column[0].split('_')
        current_chrom = first_col[1]
        
        # Write up to 1000 lines of each chromosome
        if current_chrom == temp_chrom:
            counter = counter + 1
            if counter < 1000:
                testvcf.write(line)
                
        # If a new chromosome, write a line, start counter at 0
        elif current_chrom != temp_chrom:
            counter = 0
            temp_chrom = current_chrom
            testvcf.write(line)
            
    # Include the meta lines and header line
    else: 
        testvcf.write(line)

            
testvcf.close()
vcfpath.close()