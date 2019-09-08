# import python packages
import time
import subprocess # for running an R script
import random
import math
import sys  # for getting command line arguments for chromosome # range 
import os   # for changing directory
import os.path as path  # for opening a path relative to working directory


# import classes
import SNP
import Window


# NOTE don't want b values of over about 400.
# Want calculations to be with in 5 significant figures of results from C# code
# NOTE perhaps create a user input text file where user can enter chromosome range and minimum reads number

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
    Locations = []              # Holds locations of accepted SNPs
    zraw = []                   # Holds Z scores from divergence
    #accepted_snps = []         # Number of accepted SNPs    
    outcomes = [0,0,0]          # Keeps track of [# SNPs passed, # SNPs with frequency of 1 or 0, and # SNPs carried on for B analysis]
    z_std = []                  # For standardized divergence. Use to store diverge until we get a vdiv and then go back through to calculate standardized divergence
        
    output_snps = []            # for holding accepted snps for output
    trimmed_snps = []           # for after looking at B* value
        
    min_reads = 20              # User may define option for this - seems standard but may need to be adjusted if data is nanopore sequenced... will revisit.
    Var_snp_specific = 0.0      # Used to store total variance for control and test populations
    num_snps = 0                # Number of all SNPs counter, not just accepted SNPs

    out0 = open("Results" + str(6) + ".txt", "w") # a more verbose out3
    out0.write("Chromosome\tPOS\tC AD\tControl qC_hat\tT AD\tTest qT_hat\tDivergence\n")                     
                       
    out3 = open("B" + str(6) + ".txt", "w")
    out3.write("Chromosome_Position\tB Value\tB* Value\tP Value\n")

    # STEP 2 ----------
    # Enumerate chisq file
    # TODO NOTE: probably only have to do this once in the program

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
                                              
                        # Get data from Control Sample column(s)
                        # NOTE len(info) == 5 skips over columns with no data ("./.")
                        # NOTE "j < 13" only looks at samples CA, CB, and S
                        # TODO this is a general use issue. Not all data will have 5 columns of samples
                        if len(info) == 5 and (j < 13): 
                            AD = info[1].split(",")             # AD = allele depth
                            if int(AD[0]) + int(AD[1]) > 0:    
                                c_count += 1
                                c_data.append(int(AD[0]))       # AD[0] is number of reads that support ref allele
                                c_data.append(int(AD[1]))       # AD[1] is number of reads that support alt allele
                        
                        # Get data from Test Sample column(s)
                        # NOTE "len(info) == 5" skips over columns with no data ("./.")
                        # NOTE "j > 12" only looks at samples TA and TB
                        # TODO this is a general use issue. Not all data will have 5 columns of samples.
                        if len(info) == 5 and (j > 12): 
                            AD = info[1].split(",")             # AD = allele depth
                            if int(AD[0]) + int(AD[1]) > 0:     
                                t_count += 1
                                t_data.append(int(AD[0]))       # AD[0] is number of reads that support the ref allele
                                t_data.append(int(AD[1]))       # AD[1] is number of reads that support the alt allele

                    if c_count >= 0 and t_count >= 0: 
                        qC = [0.0, 0.0]     # control group: [ref base allele depth, ref base + alt base allele depth]
                        qT = [0.0, 0.0]     # test group: [ref base allele depth, ref base + alt base allele depth]                      

                        # Control populations samples 
                        # NOTE can this be simpler? We are already getting the reference base into t_data.
                        # Why get it AGAIN to put into qC[0]? Also assign qC[1] earlier too.
                        for j in range(len(c_data) // 2):
                            m = c_data[2 * j] + c_data[2 * j + 1]
                            qC[0] += float(c_data[2 * j])    # ref base
                            qC[1] += float(m)                # ref + alt base

                        # Test populations samples
                        for j in range(len(t_data) // 2):
                            m = t_data[2 * j] + t_data[2 * j + 1]
                            qT[0] += float(t_data[2 * j])   # ref base
                            qT[1] += float(m)               # ref + alt bases

                        if (qT[1] >= min_reads and qC[1] >= min_reads):

                            # Get this first so we can tell if it actually works to take the SNP
                            qC_hat = qC[0] / qC[1]  # control groups: ref base AD / ref base + alt base AD
                            qT_hat = qT[0] / qT[1]  # test groups: ref base AD / ref base + alt base AD
                            
                            if (qC_hat != 0 and qC_hat != 1 and qT_hat != 0 and qT_hat != 1):                 
                                outcomes[2] += 1
                                Locations.append(cols[0] + "_" + cols[1])
                                num_snps += 1
                                var_C = 1.0 / qC[1]     # transformed variance of C population
                                var_T = 1.0 / qT[1]     # transformed vairance of T population

                                Var_snp_specific += (var_C + var_T)                                

                                # Calculate divergence
                                diverge = 2.0 * (math.asin(qT_hat**0.5) - math.asin(qC_hat**0.5))
                                out0.write(cols[0] + '\t' + cols[1] + '\t' + str(qC[1]) + '\t' + str(qC_hat) + '\t' + str(qT[1]) + '\t' + str(qT_hat) + '\t' + str(diverge) + '\n')
                                if random.randrange(1, 3) == 1:
                                    zraw.append(diverge)
                                else:
                                    zraw.append(-diverge)
                                    
                                z_std.append(diverge)

                                # NOTE from original code:
                                # vdiv used to be calculated here, but it was based off the var_neutral that was highly specific and precalculated in
                                # the original python script - changing the window size would have resulted in an incorrect var_neutral.
                                # The solution was to keep a seperate list of only (+) diverge and calculate it programmatically after we get the z percentiles
                                #vdiv = Var_neutral + var_C + var_T
                      
                                # We have the data we need, create a SNP and add to the output list
                                # Transform variance for each the control and test populations
                                t_C_pop = math.asin(math.sqrt(qC_hat)) # transformed variance for C population
                                t_T_pop = math.asin(math.sqrt(qT_hat)) # transformed variance for T population

                                # Create a new SNP object
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

                        # Printing file processing progress to console
                        if line_idx % 100000 == 0:
                            print(cols[0], line_idx)

    # STEP 4 ----------
    # Report on VCF file evaluation
    print("Outcomes: {}".format(outcomes))
    print("Included Number of SNPs: {}".format(sum(outcomes)))
    print("Accepted Number of SNPs: {}".format(len(output_snps)))      
    print("Sampling/genotyping variance: {}".format(Var_snp_specific / float(num_snps)))
        
    # Sort z raw list 
    ranked_z = sorted(zraw)
    n25 = ranked_z[int(num_snps / 4)]
    n50 = ranked_z[int(num_snps / 2)]
    n75 = ranked_z[int(3 * num_snps / 4)]
        
    # Report
    print("Z percentiles (without direction): {} {} {}".format(n25, n50, n75))
    print("Total variance in Z (based on IQR): {}".format(((n75 - n25) / 1.349)**2))
    print("Number skipped to do > 1 base at site: {}".format(outcomes[0]))
    print("Number of SNPs with population freq. of 1 or 0: {}".format(outcomes[1]))
    print("Number of SNPs carried for B analysis: {}".format(outcomes[2]))

    Var_neutral = ((n75 - n25) / 1.349)**2 - Var_snp_specific/float(num_snps)  # this is the bulk variance
    print("Bulk sampling and library variance {}".format(Var_neutral))

    # Run through all accepted SNPs to calculate standardized divergence ^ 2 :: B for 1 snp

    # NOTE from original code:
    # Since we only just got var neutral now we need to go back through the z_std list with only the (+) divergence and
    # finish the calculation to standardize and keep going

    for k in range(0, len(z_std)): 
        vdiv = Var_neutral + output_snps[k].c_variance + output_snps[k].t_variance
        z_std[k] = z_std[k] / (abs(vdiv)**0.5)    # NOTE do we need to add abs() here? 3/10/19

    # Storing standardized divergence
    ranked_z = sorted(z_std)
    n25 = ranked_z[int(num_snps / 4)]
    n50 = ranked_z[int(num_snps / 2)]
    n75 = ranked_z[int(3 * num_snps / 4)]
    print("Zs percentiles: {} {} {}".format(n25, n50, n75))

    print("VCF read complete, starting analysis...")



    # STEP 5: Analysis
    snploc = []     # List to hold SNP objects     
    Braw = []       # List for B values
    Bloc = []       # List to sort B

    print("Calculating B for {} SNPs...".format(outcomes[2]))

    # To grab the bp position of each of the snps in the window,
    # let's generate a new snp list to output

    all_window_SNPs = [] # a list to hold Window objects

    for k in range(window, num_snps):
        if k % window == 0 or k % window == window / 2:
                
            # Nothing is being done with this window stuff, so can skip it
            # set up window info
            new_window = Window.Window() # skip
            new_window.chromosome = output_snps[k].chromosome # skip
            new_window.window_size = k # skip

            vdiv = Var_neutral + output_snps[k].c_variance + output_snps[k].t_variance
            b = 0.0
            for j in range(window): # This is the part specifically that iterates through each window
                b += ( z_std[k-j]**2 )
                    
                new_window.snps.append(output_snps[k]) # add the snp to the list # skip
            
            all_window_SNPs.append(new_window) # TEST # skip
            output_snps[k].b_standard = b
            snploc.append(output_snps[k]) # adds the last snp of the window with value b

            Bloc.append(Locations[k])
            Braw.append(b)
            
    ranked_B = sorted(Braw)

    # Percentiles
    n25 = ranked_B[int(len(Braw) / 4)]
    n50 = ranked_B[int(len(Braw) / 2)]
    n75 = ranked_B[int(3 * len(Braw) / 4)]

    print("B Percentiles: {} {} {}".format(n25, n50, n75))
    b_skew = (n75 + n25 - 2 * n50) / (n75 - n25)
    print("B Bowley Skew: {}".format(b_skew))

    m = -1
    if b_skew > bs[0]:
        print("Too much skew")

    else:
        for j in range(1, len(bs)):
            if b_skew > bs[j]:
                m = df[j]
                jstar = j
                break

    print ("Degrees of freedom: {}".format(m))
    cIQR = percentiles[jstar][2] - percentiles[jstar][0] 
    sigB = (n75 - n25) * (2 * m)**0.5 / cIQR
    print("cIQR {} \nsigB {}".format(cIQR, sigB))

    print("Calculating B* for {} SNPs".format(len(Braw)))
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
        
    # TODO write the test data from snp window -- ASK WHAT THIS MEANS
    # Allows to keep track of whole window, otherwise losing 66 % of window details
       
    # Remove all where there is no b* i.e. b_star is less than or equal to 0      
    finalSNPlist = [snp for snp in snploc if snp.b_star > 0]
      
    # Close files
    out3.close()
    out0.close()

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
    
    print("\n{} SNPs removed below FDR threshold leaving {}".format(pre_removal_length, len(sorted_list)))

    # Now show the sig b*
    min_SNP = min(sorted_list, key = lambda snp: snp.b_star)
    print("\nSignificant B* (the min B* after removing SNPs over threshold) {}".format(min_SNP.b_star)) # TODO this prints out an address, needs to do a string

    return sorted_list # return a list of SNPs

# ------------------------------------ END FDR ANALYSIS ------------------------------------------------ #

# ------------------------------------ BEGIN MAIN ------------------------------------------------------ #


# STEP 1 ----------
#  - open the vcf and chisq filepaths
#  - NOTE for testing, enter the correct paths
start_recording = time.time()
try:
    #vcf_path = open(path.abspath(path.join(sys.path[0], "test_files/REDUCED_ali.vcf")), "r")
    vcf_path = open(path.abspath(path.join(sys.path[0], "test_files/Ali_w_767.vcf")), "r")
    chisq_path = open(path.abspath(path.join(sys.path[0], "support_files/chisq.txt")), "r")

except IOError:
    print("Error opening file paths")

# STEP 2 ----------
#  - Get chromosome range for reading vcf file
if len(sys.argv) == 3:
    lowerlimit = int(sys.argv[1]) # "1"  for testing
    upperlimit = int(sys.argv[2]) # "14" for testing
else:
    print("Invalid arguments for chromosome range. Please put beginning and ending chromosome numbers.")


# Change directory to results to save resulting files in
os.chdir("results")

#   STEP 3 ----------
#  - run the vcf parser to determine b value for a SNP window of 1 (s = 1) 
#  - creates a snp list file called B1_new.txt, which genwin will read later
start_time = time.time()
print("Starting B processing at: {}".format(time.asctime(time.localtime())))
    
write_results = open("B1_new.txt", "w+")
write_results.write("CHR\tBP\tB\n")
    
snp_list = VCF_Analyzer(1, vcf_path, chisq_path, lowerlimit, upperlimit)
    
for snp in snp_list:
    write_results.write("" + str(snp.chromosome) + '\t' + str(snp.basepair) + '\t' + str(snp.b_standard) + "\n") 
elapsed_time = time.time() - start_time
print("B processing time: {} seconds".format(elapsed_time))

# Close files
write_results.close()
vcf_path.close()
chisq_path.close()
    
# STEP 4 ----------
#  - calls the R script, GenWin 
#  - this reads the B1_new.txt file that was created above
#  - GenWin creates a file called "splinewindows.txt"
start_time = time.time()
print("Beginning R execution at: {}".format(time.asctime(time.localtime())))
    
# NOTE Potential Issue with calling R:
#  - Having a hard coded path to R may be a problem for users who have installed their R in a different directory
#  - if user has not set a mirror for their R, then an error will appear that the user is trying to install packages without setting a mirror

rcommand_path = "C:/Program Files/R/R-3.5.1/bin/Rscript.exe" 
    
rscript_path = path.abspath(path.join(sys.path[0], "support_files/GenWin_script_12_29_2016.R"))

# AT HOME
#rscript_path = "C:\Users\Heathro\OneDrive\Mimulus\mg-gap-py\mg-gap\support_files\GenWin_script_12_29_2016.R" 

# AT CWU
#rscript_path = "C:/Users/gammonh/OneDrive/Mimulus/mg-gap-py/mg-gap/support_files/GenWin_script_12_29_2016.R"
subprocess.call([rcommand_path, rscript_path])
elapsed_time = time.time() - start_time
print("R exited successfully.\nRun time: {} seconds".format(elapsed_time))
    
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
splinepath = path.abspath(path.join(sys.path[0], "results\splinewindows.txt")) 

# AT HOME
#splinepath = "C:\Users\Heathro\OneDrive\Mimulus\mg-gap-py\mg-gap\splinewindows.txt" 

# AT CWU
#splinepath = "C:/Users/gammonh/OneDrive/Mimulus/mg-gap-py/mg-gap/splinewindows.txt"
print("GenWin file found, obtaining median window value...")
windows = []    # List to hold window sizes to find median value
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
		print("The median widow size is: {}".format(median))
except:
	print("Error opening splinewindows.txt")



# STEP 6 ----------
#  - reopen vcf and chisq paths
#  - feed the median value back through b processing, then b* 
#  - need to hold on to the data for FDR 
#  - make a tab .csv
#median = 6 # TODO remove this - for testing only
try:
    #vcf_path = open(path.abspath(path.join(sys.path[0], "test_files/REDUCED_ali.vcf")), "r")
    vcf_path = open(path.abspath(path.join(sys.path[0], "test_files/Ali_w_767.vcf")), "r")
    chisq_path = open(path.abspath(path.join(sys.path[0], "support_files/chisq.txt")), "r")

except IOError:
	print("Error opening file paths")

if median > 0:
	snpList = VCF_Analyzer(median, vcf_path, chisq_path, lowerlimit, upperlimit)
	print("Re-analyzing for B* based on median window size {} @ {}".format(median, time.asctime(time.localtime())))

	start_time = time.time()
	print("Starting FDR sorting processing at: {}".format(time.asctime(time.localtime())))

	# Start the FDR process
	# This should be a user config variable in the future if they would prefer any different settings
	# TODO maybe create user input for fdr_input value
	fdr_input = 0.05
	print("Running FDR analysis at {}...".format(fdr_input))
	fdrlist = process(snpList, fdr_input)

	# Save this fdrlist to a "tab" separated values, use a .csv, file with headers 
	# TODO this is the same code as in step 2 - good candidate for a write snp list method
	write_results = open("FDR_analysis.txt", "w+")
	write_results.write("CHR\tBP\tBstar\tRawP\n")
	for snp in snpList:
		write_results.write("" + str(snp.chromosome) + '\t' + str(snp.basepair) + '\t' + str(snp.b_star) + "\t" + str(snp.raw_p) + "\n")
	elapsed_time = time.time() - start_time
	print("Sort FDR processing time: {} seconds".format(elapsed_time))

	# Start getting the annotations !! Can stop here for the purposes of this program. 
	# TODO Something that may be of interest later. Connect to pythosome site to get Gene, description, and RNA SEQ data.

# Close files
vcf_path.close()
chisq_path.close()
write_results.close()

# END WHOLE PROGRAM RUN TIME
end_recording = time.time()
whole_recording = end_recording - start_recording
print("Total program run time: {}".format(whole_recording))

# Write whole program run time to file
#time_record = open("py_time_record.txt", "a")
#time_record.write("\n{}".format(whole_recording))
#time_record.close()
# ------------------------------------ END MAIN ------------------------------------------------------ #