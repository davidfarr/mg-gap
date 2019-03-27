# import python packages
import time
import subprocess # for running an R script
import random
import math
import sys # for getting command line arguments for chromosome # range 


# import classes
import SNP
import Window


# NOTE don't want b values of over about 400.
# Want calculations to be with in 5 significant figures.
# NOTE perhaps create a user input text file where user can enter chromosom range and minimum reads number

    # ------------------------------------- BEGIN VCF ANALYZER --------------------------------------------- #
  
def VCF_Analyzer(window, vcfpath, chisq_path, lowerlimit, upperlimit):
    """
    Function: VCF_Analyzer

    Description: Reads the VCF. Create a list of SNP's with their info
    in a method and return the "list" of SNPs back to the main program file.

    @param window (integer): the median window value to look at SNPs
    @param vcfpath (string): file path to vcf file
    @param chisq_path (string): file path to chi squared file
    @param lowerlimit (int): lower limit of chromosome analysis range
    @param upperlimit (int): upper limit of chromosome analysis range
    @return finalSNPlist (SNP list): a list of SNP objects with b* values

    DEVELOPMENT NOTES: 
     - it is very long. Maybe we can break it up into separate, smaller methods.
     - when translating, replaced S with window, src with vcfpath, in1 with chisq_path
     - maybe move some of the STEPS into their own functions for ease of reading and testing
     - a couple counters, total_counter and skip_counter are used in the C#, do we want to put this
     back in?
     - there's a note at the end of this code too, but do we want to write the test data from snp window?
    """
    # STEP 1 ----------
    # Set up counters, lists, and variables
    Locations = []              # ?
    zraw = []                   # Holds Z scores from divergence
    accepted_snps = []          # Number of accepted SNPs    
    outcomes = [0,0,0]          # Keeps track of [# SNPs passed, # SNPs with frequency of 1 or 0, and # SNPs carried on for B analysis]
    z_std = []                  # For standardized divergence. Use to store diverge until we get a vdiv and then go back through to calculate standardized divergence
        
    output_snps = []            # from C# code, for holding output snps
    trimmed_snps = []           # from C# code, for after B*
        
    min_reads = 20              # User may define option for this - seems standard but may need to be adjusted if data is nanopore sequenced... will revisit.
    Var_snp_specific = 0.0      # ?
    num_snps = 0                # Number of all SNPs counter, not just accepted SNPs

    out0 = open("Results" + str(6) + ".txt", "w") # a more verbose out3
    out3 = open("B" + str(6) + ".txt", "w")

    # STEP 2 ----------
    # Enumerate chisq file
    print("Assembling chi square apparatus...")
    bs = [] 
    df = []
    percentiles=[[] for i in range(100)]
    for line_idx, line in enumerate(chisq_path):
        cols = line.replace('\n', '').split('\t')    
        df.append(float(cols[0]))
        for j in range(9):
            percentiles[line_idx].append(float(cols[j + 1]))
        rx = (percentiles[line_idx][2] + percentiles[line_idx][0] - 2 * percentiles[line_idx][1]) / (percentiles[line_idx][2] - percentiles[line_idx][0])
        bs.append(rx)

    # STEP 3 ----------
    # Evaluate contents of each line of the VCF input file          
    for line_idx, line in enumerate(vcfpath):
        cols = line.replace('\n', '').split('\t') # Get rid of newlines at the end of each line, then split each line on the tab 
        if len(cols) < 2: 
            continue # Skip the meta section (file description, definitions, contig info, etc)
        elif cols[0] == "#CHROM":       # this is the header row for the variant sample data
            for i in range(len(cols)):
                print(i, cols[i])

        else:              
            # Get the actual chromosome number out of the first column by removing the "sNNafold_"
            # NOTE general use issue, not all VCF files will have the "sNNffold_"
            # TODO fix this ^ - read the number only from this column, ignore non digit characters
            chromosome = cols[0]        # sNNffold_x, where x is the number
            chromosome = chromosome[9:] # x (we got rid of everything so it's just the number)

            # Analyze the chromosomes in the user-defined range
            if lowerlimit <= int(chromosome) and int(chromosome) <= upperlimit:   
                position = int(cols[1]) # grab the base pair position
                ref_base = cols[3]      # the reference base allele
                alt_base = cols[4]      # the alternate base allele

                # Check for multiple bases
                if len(alt_base) > 1:   # Too many bases at site
                    outcomes[0] += 1    # Pass this SNP
                else:
                    # Set up more counters
                    c_count = 0      # control sample count
                    c_data = []      # control Allele Depth data
                    t_count = 0      # test sample count
                    t_data = []      # test Allele Depth data
                                                        
                    # NOTE are skipping a sample column, not all reserachers will want to skip this column
                    # TODO re-create ali vcf without the 767, make note in user manual that users must create vcf's that 
                    # contain only samples they want to look at
                    for j in range(10, len(cols)):      # ali_w_767.vcf... keep in mind that the data *structure* difference between ali.vcf and w/767 is the 767 bam file is the first bam column index, so must add +1 to all the indexes that referenced bam file cols. Skip 767 which is column 10.
                        info = cols[j].split(':')       # split the info column that has AD on col[1]
                                              
                        # len(info) == 5 skips over columns with no data ("./."), j < 13 only looks at samples CA, CB, and S
                        # TODO this is a general use issue. Not all data will have 5 columns of samples
                        if len(info) == 5 and (j < 13): 
                            AD = info[1].split(",")             # AD = allele depth
                            if int(AD[0]) + int(AD[1]) > 0:     # AD[0] is number of reads that support ref allele, AD[1] is number of reads that support alt allele
                                c_count += 1
                                c_data.append(int(AD[0]))
                                c_data.append(int(AD[1]))
                        
                        # len(info) == 5 skips over columns with no data ("./."), j > 12 only looks at samples TA and TB
                        # TODO this is a general use issue. Not all data will have 5 columns of samples
                        if len(info) == 5 and (j > 12): 
                            AD = info[1].split(",")             # AD = allele depth
                            if int(AD[0]) + int(AD[1]) > 0:     # AD[0] is number of reads that support the ref allele, AD[1] is number of reads that support the alt allele
                                t_count += 1
                                t_data.append(int(AD[0]))
                                t_data.append(int(AD[1]))

                    if c_count >= 0 and t_count >= 0: 
                        qC = [0.0, 0.0]     # control group: ref base allele depth, ref+alt base allele depth
                        qT = [0.0, 0.0]     # test group: ref base allele depth, ref+alt base allele depth                        

                        # Control Samples 
                        # NOTE can this be simpler? We are already getting the reference base into t_data, why get it AGAIN to put into qC[0]? Also assign qC[1] earlier too
                        for j in range(len(c_data) // 2):
                            m = c_data[2 * j] + c_data[2 * j + 1]
                            qC[0] += float(c_data[2 * j])    # ref base
                            qC[1] += float(m)                # ref + alt base

                        # Test Samples
                        for j in range(len(t_data) // 2):
                            m = t_data[2 * j] + t_data[2 * j + 1]
                            qT[0] += float(t_data[2 * j])   # ref base
                            qT[1] += float(m)               # ref + alt bases

                        if (qT[1] >= min_reads and qC[1] >= min_reads):

                            # Get this first so we can tell if it actually works to take the SNP
                            qC_hat = qC[0] / qC[1]  # control group: ref base AD / ref+alt base AD
                            qT_hat = qT[0] / qT[1]  # test group: ref base AD / ref+alt base AD
                            
                            if (qC_hat != 0 and qC_hat != 1 and qT_hat != 0 and qT_hat != 1):
                                # columns: scaffold, position, ref allele count 1, total count 1,ref allele count 2, total count 2
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

    # STEP 4 ----------
    # Report on VCF file evaluation
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
        z_std[k] = z_std[k] / (abs(vdiv)**0.5)    # NOTE do we need to add abs() here? 3/10/19

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
        print("Too much skew")

    else:
        for j in range(1, len(bs)):
            if b_skew > bs[j]:
                m = df[j]
                jstar = j
                break

    print ("Degrees of freedom: ", m)
    cIQR = percentiles[jstar][2] - percentiles[jstar][0] 
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
        
    # TODO write the test data from snp window
    # Allows to keep track of whole window, otherwise losing 66 % of window details
       
    # Remove all where there is no b* i.e. b_star is less than or equal to 0      
    finalSNPlist = [snp for snp in snploc if snp.b_star > 0]
      
    print("VCF analysis complete.")
    return finalSNPlist # return a list of SNPs


# ------------------------------------- END VCF ANALYZER ----------------------------------------------- #


# ------------------------------------- BEGIN FDR ANALYSIS --------------------------------------------- #

def process(bs_list, fdr_selected):
    """
    Function: FDR

    Description: 
    The FDR process is:
        1. Sort B* list by p-value
        2. Add new property FDR where FDR = 0.1 * index of the ranked SNP / (count of SNPs in the file / 2)
        3*. FDR can be changed... 0.1 above is FDR of 10 and 0.05 is FDR 5

    @param bs_list (SNP object list): the list of SNPs with b* values
    @param fdr_selected - TODO need to learn more about FDR
    @return sorted_list (SNP object list): a list of SNPs with significant b* values

    Editing thoughts for later: TODO request to use FDR analyzer in R instead.
    """  
    sorted_list = sorted(bs_list, key = lambda snp: snp.raw_p)
        
    rank_assignment = 1
    for snp in sorted_list:
        snp.fdr_rank = rank_assignment
        rank_assignment += 1
         
    num_window = rank_assignment / 2
    for snp in sorted_list:
        snp.threshhold_value = (fdr_selected * snp.fdr_rank) / num_window 

    bs_capacity = len(bs_list)
    for i in range(0, len(bs_list)):
        bs_list[i].threshold_value = fdr_selected * (i + 1) / (bs_capacity / 2)

    # TEST
    # remove SNPs that have a raw_p value > threshold_value   
    pre_removal_length = len(sorted_list)
    sorted_list = [snp for snp in sorted_list if snp.raw_p <= snp.threshold_value]
    
    print("\n%s SNPs removed below FDR threshold leaving %s" % (pre_removal_length, len(sorted_list)) )

    # Now show the sig b*
    print("\nSignificant B* (the min B* after removing SNPs over threshold) %s" % min(sorted_list, key = lambda snp: snp.b_star)) # TODO this prints out an address, needs to do a string

    return sorted_list # return a list of SNPs

# ------------------------------------ END FDR ANALYSIS ------------------------------------------------ #

# ------------------------------------ BEGIN MAIN ------------------------------------------------------ #

# STEP 1 ----------
#  - open the vcf and chisq filepaths
#  - NOTE for testing, enter the correct paths
try:
    # vcf_path = open("D:/MG_GAP/mg-gap/mg-gap/mg-gap-py/mg-gap/test_files/REDUCED_ali.vcf", "rU") TODO remove U, deprecated 
	vcf_path = open("C:/Users/gammonh/Documents/GitHub/mg-gap/mg-gap/mg-gap-py/mg-gap/test_files/REDUCED_ali.vcf", "r") # Grad lab office PC
    # chisq_path = open("D:/MG_GAP/mg-gap/mg-gap/mg-gap-py/mg-gap/support_files/chisq.txt", "rU") TODO remove U, deprecated
	chisq_path = open("C:/Users/gammonh/Documents/GitHub/mg-gap/mg-gap/mg-gap-py/mg-gap/support_files/chisq.txt", "r") # Grad lab office PC
except:
    print("Error opening file paths")

# STEP 2 ----------
#  - Get chromosome range for reading vcf file
if len(sys.argv) == 3:
	lowerlimit = int(sys.argv[1])
	upperlimit = int(sys.argv[2])
else:
    print("Invalid arguments for chromosome range. Please put beginning and ending chromosome numbers.")

# STEP 3 ----------
#  - run the vcf parser to determine b value for a SNP window of 1 (s = 1) 
#  - creates a snp list file called B1_new.txt, which genwin will read later
start_time = time.time()
print("Starting B processing at :", start_time)

write_results = open("B1_new.txt", "w+")
write_results.write("CHR\tBP\tB\n")
snp_list = VCF_Analyzer(1, vcf_path, chisq_path, lowerlimit, upperlimit)
for snp in snp_list:
    write_results.write("" + str(snp.chromosome) + '\t' + str(snp.basepair) + '\t' + str(snp.b_standard) + "\n") 
elapsed_time = time.time() - start_time
print("B processing time: ", elapsed_time)
write_results.close()
vcf_path.close()
chisq_path.close()

# STEP 4 ----------
#  - calls the R script, GenWin 
#  - this reads the B1_new.txt file that was created above
#  - GenWin creates a file called "splinewindows.txt"
start_time = time.time()
print("Beginning R execution at", start_time)

# NOTE Potential Issue with calling R:
#  - Having a hard coded path to R may be a problem for users who have installed their R in a different directory
#  - if user has not set a mirror for their R, then an error will appear that the user is trying to install packages without setting a mirror
rcommand_path = "C:/Program Files/R/R-3.5.1/bin/Rscript.exe" 
#rscript_path = "D:/MG_GAP/mg-gap/mg-gap/mg-gap-py/mg-gap/support_files/GenWin_script_12_29_2016.R" # Personal home PC
rscript_path = "C:/Users/gammonh/Documents/GitHub/mg-gap/mg-gap/mg-gap-py/mg-gap/support_files/GenWin_script_12_29_2016.R" # Grad lab office PC
subprocess.call([rcommand_path, rscript_path])
elapsed_time = time.time() - start_time
print("R exited successfully.\nRun time: ", elapsed_time)
 

# STEP 5 ----------
#  - calculate the median from the "splinewindows.txt" file created by GenWin
#  - run the b processing at that window and then b* processing 
#  - check if we recognize that the file has now been put there
"""
Guide to splinewindows format
   * 0 = CHRcol (snnaffold #)
   * 1 = Window start bp position
   * 2 = window end bp position
   * 3 = SNP count (number of snps in the window)
   * 4 = Mean Y
   * 5 = W statistic
   * first row is headers, skip
"""
#splinepath = "D:/MG_GAP/mg-gap/mg-gap/mg-gap-py/mg-gap/splinewindows.txt" # Personal home PC
splinepath = "C:/Users/gammonh/Documents/GitHub/mg-gap/mg-gap/mg-gap-py/mg-gap/splinewindows.txt" # Grad lab office PC
print("GenWin file found, obtaining median window value...")
windows = []
try:
    with open(splinepath, "r") as sp:
        sp.readline() # read past the header line
        for line in sp:
            if ("CHRcol" not in line):
                cols = line.replace("\n", "").split("\t")
                windows.append(int(cols[3])) 
    if len(windows) == 0:
        print("No window sizes to calculate!")
    else:
        windows.sort()
        middle = len(windows) // 2
        if (len(windows) % 2 != 0) :
            median = windows[middle]
        else:
            median = (windows[middle] + windows[middle - 1]) // 2
        print("The median widow size is: ", median)
except:
    print("Error opening splinewindows.txt")


# STEP 6 ----------
#  - reopen vcf and chisq paths
#  - feed the median value back through b processing, then b* 
#  - need to hold on to the data for FDR 
#  - make a tab .csv
#median = 6 # TODO remove this - for testing only
try:
    # vcf_path = open("D:/MG_GAP/mg-gap/mg-gap/mg-gap-py/mg-gap/test_files/REDUCED_ali.vcf", "rU") TODO remove U, deprecated
	vcf_path = open("C:/Users/gammonh/Documents/GitHub/mg-gap/mg-gap/mg-gap-py/mg-gap/test_files/REDUCED_ali.vcf", "r") # Grad lab office PC
    # chisq_path = open("D:/MG_GAP/mg-gap/mg-gap/mg-gap-py/mg-gap/support_files/chisq.txt", "rU") TODO remove U, deprecated
	chisq_path = open("C:/Users/gammonh/Documents/GitHub/mg-gap/mg-gap/mg-gap-py/mg-gap/support_files/chisq.txt", "r") # Grad lab office PC
except:
    print("Error opening file paths")

if median > 0:
    snpList = VCF_Analyzer(median, vcf_path, chisq_path, lowerlimit, upperlimit)
    print("Re-analyzing for B* based on median window size ", median, " @ ", time.ctime)

    start_time = time.time()
    print("Starting FDR sorting processing at :", start_time)

    # Start the FDR process
    # This should be a user config variable in the future if they would prefer any different settings
    # TODO maybe create user input for fdr_input value
    fdr_input = 0.05
    print("Running FDR analysis at ", fdr_input, "...")
    
    fdrlist = process(snpList, fdr_input)

    # Save this fdrlist to a "tab" separated values, use a .csv, file with headers 
    # TODO this is the same code as in step 2 - good candidate for a write snp list method
    write_results = open("FDR_analysis.txt", "w+")
    write_results.write("CHR\tBP\tB\n")
    for snp in fdrlist:
        write_results.write("" + str(snp.chromosome) + '\t' + str(snp.basepair) + '\t' + str(snp.b_standard) + "\n")
    elapsed_time = time.time() - start_time
    print("Sort FDR processing time: ", elapsed_time)

  
    # Start getting the annotations !! Can stop here for the purposes of 
    # this program. 

    # TODO Something that may be of interest later. Connect to pythosome site to get Gene, description, and RNA SEQ data.


    # ------------------------------------ END MAIN ------------------------------------------------------ #