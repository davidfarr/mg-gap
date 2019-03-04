"""
Class VCF_Analyzer

Description: Reads the VCF. Create a list of SNP's with their info
in a method and return the "list" of SNPs back to the main program file.

NOTES: 
 - this class contains only one method. No fields exist for this class 
 when it is instantiated. Perhaps move method to the main program.py 
 file as a method outside of main? Also, it is very long. Maybe we 
 can break it up into separate, smaller methods.
 - when translating, replaced S with window, src with vcfpath, in1 with chisq_path
 - maybe move some of the STEPS into their own functions for ease of reading and testing
 - a couple counters, total_counter and skip_counter are used in the C#, do we want to put this
 back in?
 - there's a note at the end of this code too, but do we want to write the test data from snp window?
"""
# class imports
import SNP
import Window
import math
import random

#class VCF_Analyzer:

   
def SNP_list(window, vcfpath, chisq_path): # need to add keyword "self" if adding actual class properties
        
    # STEP 1: Set up counters, lists, and variables
    LineID = []                 # For holding column names (headers) from VCF file
    Locations = []              # ?
    zraw = []                   # Holds Z scores from divergence
    accepted_snps = []          # Number of accepted SNPs    
    outcomes = [0,0,0]          # ?
    z_std = []                  # For standardized divergence. Use to store diverge until we get a vdiv and then go back through to calculate standardized divergence
        
    output_snps = []            # from C# code, for holding output snps
    trimmed_snps = []           # from C# code, for after B*
        
    min_reads = 20              # User may define option for this - seems standard but may need to be adjusted if data is nanopore sequenced... will revisit.
    depth_cc = 0                # ?
    raw_read_count = 0          # ?
    Var_snp_specific = 0.0      # ?
    num_snps = 0                # Number of all SNPs counter, not just accepted SNPs

    outkk = open("outkkTEST.yut", "w") # left this in, but no current use in application TODO ask about this
    out0 = open("Results" + str(6) + ".txt", "w") #a more verbose out3
    out3 = open("B" + str(6) + ".txt", "w")

    # STEP 2: Enumerate chisq file
    print("Assembling chi square apparatus...")
    bs = [] 
    df = []
    percentiles=[[] for i in range(100)]
    for line_idx, line in enumerate(chisq_path):
        cols = line.replace('\n', '').split('\t')    
        df.append(float(cols[0]))
        # bs.append(float(cols[1]))
        for j in range(9):
            percentiles[line_idx].append(float(cols[j+1]))
        rx = (percentiles[line_idx][2] + percentiles[line_idx][0] - 2 * percentiles[line_idx][1]) / (percentiles[line_idx][2] - percentiles[line_idx][0])
        bs.append(rx)

    # STEP 3: Evaluate contents of each line of the VCF input file            
    for line_idx, line in enumerate(vcfpath):
        cols = line.replace('\n', '').split('\t') # Get rid of newlines and empty spaces, then use the tab delimiter for array
        # Should be contig row... these list contigs and have no data for use
        if len(cols) < 2:  # Sorts out any 2-col contigs
            continue # Skip the header row
        # This is indicative of the first row of the actual data set. Print out every bam file name
        elif cols[0] == "#CHROM":       # this is the header row (but not the first row)
            for i in range(len(cols)):
                print(i,cols[i])
                if i > 8:
                    LineID.append(cols[i])
        else: # this is if len(cols) >= 2, NOTE in C# code this is just length > 2
               
            # Get the actual chromosome number out of the first column
            chromosome = cols[0]        # snnafold_x
            chromosome = chromosome[9:] # x (we got rid of everything so it's just the number)
                
            if int(chromosome) < 15:    # Many more contigs than chromosomes - we're only looking between chr 1 and 14 for meaningful data.
                scaff = cols[0].split('_')
                position = int(cols[1]) # grab the base pair
                ref_base = cols[3]      # the reference base allele
                alt_base = cols[4]      # the alternate base allele

                # Check for multiple bases
                if len(alt_base) > 1:   # Too many bases at site
                    outcomes[0] += 1    # Pass this SNP
                else:
                    # Set up more counters
                    C_count = 0
                    C_dat = []
                    T_count = 0
                    T_dat = []

                    for j in range(10, len(cols)):      # ali_w_767.vcf... keep in mind that the data *structure* difference between ali.vcf and w/767 is the 767 bam file is the first bam column index, so must add +1 to all the indexes that referenced bam file cols. Skip 767 which is column 10.
                        info = cols[j].split(':')       # split the info column that has AD on col[1]
                        if len(info) == 5 and (j < 13): # ali_w_767.vcf
                            AD = info[1].split(",")
                            if int(AD[0]) + int(AD[1]) > 0:
                                C_count += 1
                                C_dat.append(int(AD[0]))
                                C_dat.append(int(AD[1]))
                        if len(info) == 5 and (j > 12): # ali_w_767.vcf
                            AD = info[1].split(",")
                            if int(AD[0]) + int(AD[1]) > 0:
                                T_count += 1
                                T_dat.append(int(AD[0]))
                                T_dat.append(int(AD[1]))

                    if C_count >= 0 and T_count >= 0:
                        qC = [0.0, 0.0]
                        qT = [0.0, 0.0]                              

                        for j in range(len(C_dat) // 2):
                            m = C_dat[2 * j] + C_dat[2 * j + 1]
                            qC[0] += float(C_dat[2 * j])    # ref base
                            qC[1] += float(m)               # ref + alt base
                        for j in range(len(T_dat) // 2):
                            m = T_dat[2 * j] + T_dat[2 * j + 1]
                            qT[0] += float(T_dat[2 * j])
                            qT[1] += float(m)

                        if (qT[1] >= min_reads and qC[1] >= min_reads):
                            # Get this first so we can tell if it actually works to take the SNP
                            qC_hat = qC[0] / qC[1]
                            qT_hat = qT[0] / qT[1]
                            if (qC_hat != 0 and qC_hat != 1 and qT_hat != 0 and qT_hat != 1):
                                # columns: scaffold, position, ref allele count 1, total count 1,ref allele count 2, total count 2
                                outkk.write(cols[0] + "\t"+cols[1] + "\t" + str(qC[0]) + "\t"+str(qC[1]) + "\t"+str(qT[0]) + "\t" + str(qT[1]) + "\n")
                                outcomes[2] += 1
                                Locations.append(cols[0] + "_" + cols[1])
                                num_snps += 1
                                var_C = 1.0 / qC[1]     # transformed variance of C
                                var_T = 1.0 / qT[1]     # for T

                                Var_snp_specific += (var_C + var_T)
                                #vdiv = Var_neutral + var_C + var_T

                                # Calculate divergence, NOTE do we want to add zranked here?
                                diverge = 2.0 * (math.asin(qT_hat**0.5) - math.asin(qC_hat**0.5))
                                out0.write(cols[0] + '\t' + cols[1] + '\t' + str(qC[1]) + '\t' + str(qC_hat) + '\t' + str(qT[1]) + '\t' + str(qT_hat) + '\t' + str(diverge) + '\n')
                                if random.randrange(1, 3) == 1:
                                    zraw.append(diverge)
                                else:
                                    zraw.append(-diverge)
                                    
                                z_std.append(diverge)

                                #Major change:
                                #vdiv used to be calculated here, but it was based off the var_neutral that was highly specific and precalculated in
                                #the original python script - changing the window size would have resulted in an incorrect var_neutral.
                                #The solution was to keep a seperate list of only (+) diverge and calculate it programmatically after we get the z percentiles

                                accepted_snps.append([cols[0] + "_" + cols[1], var_C, var_T])
                                    
                                # NOTE this block was added following the C# code                        
                                # We have the data we need, create a SNP and add to the output list
                                t_C_pop = math.asin(math.sqrt(qC_hat)) # transformed variance for C population
                                t_T_pop = math.asin(math.sqrt(qT_hat)) # transformed variance for T population

                                acceptedSNP = SNP.SNP()
                                acceptedSNP.chromosome = chromosome
                                acceptedSNP.basepair = position
                                acceptedSNP.c_variance = var_C
                                acceptedSNP.t_variance = var_T
                                acceptedSNP.transformed_c_variance = t_C_pop 
                                acceptedSNP.transformed_t_variance = t_T_pop 
                                acceptedSNP.old_identifier = cols[0]
                                output_snps.append(acceptedSNP)                                   

                            else:
                                outcomes[1] += 1

                        if line_idx % 100000 == 0:
                            print(cols[0], line_idx)

    # STEP 4: Report on VCF file evaluation
    print("Outcomes: ", outcomes)
    print("Included: Number of SNPs ", num_snps)
    print("Accepted: Number of SNPs ", len(output_snps)) # from C# code        
    print("Sampling/genotyping variance: ", Var_snp_specific / float(num_snps))
        
    # Sort z list (only the one actually for ranking)
    ranked_z = sorted(zraw)
    n25 = ranked_z[int(num_snps / 4)]
    n50 = ranked_z[int(num_snps / 2)]
    n75 = ranked_z[int(3 * num_snps / 4)]
        
    # Report
    print("Z percentiles (without direction): ", n25, n50, n75)
    print("Total variance in Z (based on IQR): ",((n75 - n25) / 1.349)**2)
    print("Number skipped to do > 1 base at site: %s\n"          % outcomes[0] +
            "Number of SNPs with population freq. of 1 or 0: %s\n" % outcomes[1] +
            "Number of SNPs carried for B analysis: %s"            % outcomes[2] )
    Var_neutral = ((n75 - n25) / 1.349)**2 - Var_snp_specific/float(num_snps)  # this is the bulk variance
    print("Bulk sampling and library variance ", Var_neutral)

    # Run through all accepted SNPs to calculate standardized divergence ^ 2 :: B for 1 snp

    # Major Change:
    # Since we only just got var neutral now we need to go back through the z_std list with only the (+) divergence and
    # finish the calculation to standardize and keep going

    for k in range(0, len(z_std)):
        vdiv = Var_neutral + accepted_snps[k][1] + accepted_snps[k][2]
        z_std[k] = z_std[k] / (vdiv**0.5)

    ranked_z = sorted(z_std)
    n25 = ranked_z[int(num_snps / 4)]
    n50 = ranked_z[int(num_snps / 2)]
    n75 = ranked_z[int(3 * num_snps / 4)]
    print("Zs percentiles ", n25, n50, n75)

    print("VCF read complete, starting analysis...")


    # STEP 5: Analysis
    snploc = [] # a list to hold SNP objects

    # List for B and a list to sort B
    Braw = []
    Bloc = []

    print("Calculating B for %s SNPs..." % outcomes[2])

    # To grab the bp position of each of the snps in the window,
    # let's generate a new snp list to output

    all_window_SNPs = [] # a list to hold Window objects

    for k in range(window, num_snps):
        if k % window == 0 or k % window == window / 2:
                
            # set up window info (added in from C# code)
            new_window = Window.Window()
            new_window.chromosome = output_snps[k].chromosome
            new_window.window_size = k

            vdiv = Var_neutral + accepted_snps[k][1] + accepted_snps[k][2]
            b = 0.0
            for j in range(window): # This is the part specifically that iterates through each window
                b += ( z_std[k-j]**2 )
                    
                new_window.snps.append(output_snps[k]) # add the snp to the list
            
            # from C# code
            all_window_SNPs.append(new_window) # TEST
            output_snps[k].b_standard = b
            snploc.append(output_snps[k]) # adds the last snp of the window with value b

            Bloc.append(Locations[k])
            Braw.append(b)
            
    ranked_B = sorted(Braw)

    # Percentiles
    n25 = ranked_B[int(len(Braw) / 4)]
    n50 = ranked_B[int(len(Braw) / 2)]
    n75 = ranked_B[int(3 * len(Braw) / 4)]

    print("B Percentiles: ", n25, n50, n75)
    b_skew = (n75 + n25 - 2 * n50) / (n75 - n25)
    print("B Bowley Skew: ", b_skew)

    m = -1
    if b_skew > bs[0]:
        print("Too much skew") # in this if statement, jstar never gets initialized, this is a problem for the program to continue, need to add error checking here

    else:
        for j in range(1, len(bs)):
            if b_skew > bs[j]:
                m = df[j]
                jstar = j
                break

    print ("Degrees of freedom: ", m)
    cIQR = percentiles[jstar][2] - percentiles[jstar][0] # TODO error here
    sigB = (n75 - n25) * (2 * m)**0.5 / cIQR
    print("cIQR %s \nsigB %s" % (cIQR, sigB))

    print("Calculating B* for %s SNPs" % len(Braw))
    for j in range(len(ranked_B)):
        bx = Braw[j]
        bs = m + (bx - window) * ((2 * m)**0.5) / sigB        

        if bs < percentiles[jstar][3]: # p greater than 0.05
            p = 0.5
        elif bs > percentiles[jstar][8]:
            p = 5.0 * 10**(1.0 - 8.5) # p less than table min
        else:
            for k in range(4,9):
                if bs > percentiles[jstar][k - 1] and bs <= percentiles[jstar][k]:
                    dx = (bs - percentiles[jstar][k - 1])/(percentiles[jstar][k] - percentiles[jstar][k - 1])
                    x = k - 1 + dx
                    p = 5.0 * 10**(1.0 - x)
                    break
        if j > 0:
            midpoint = str(Bloc[j - 1])
            # from C# code
            snploc[j - 1].b_star = bs
            snploc[j - 1].raw_p = p;
            snploc[j - 1].b_standard = bx;

        else:
            midpoint = str(Bloc[j])
            # from C# code
            snploc[j].b_star = bs
            snploc[j].raw_p = p;
            snploc[j].b_standard = bx;

        out3.write(midpoint + '\t' + str(bx) + '\t' + str(bs) + '\t' + str(p) + '\n')
        
    # TODO ask if we want to write the test data from snp window thingymabob?
    # Allows to keep track of whole window, otherwise losing 66 % of window details
    # Yes do this, see C# code. 
       
    # Remove all where there is no b* i.e. b_star is less than or equal to 0        
    #finalSNPlist = list(filter((snploc.b_star <= 0).__ne__, snploc))     
     
    finalSNPlist = [snp for snp in snploc if snp.b_star > 0]
      
    print("VCF analysis complete.")
    return finalSNPlist # return a list of SNPs