# import python packages
import time
import subprocess # for running an R script

# import classes
import VCF_Analyzer
import FDR



# STEP 1 ----------
#  - open the vcf and chisq filepaths
#  - NOTE for testing, enter the correct paths
try:
    vcf_path = open("C:\Users\Public\Ali_w_767.vcf", "rU")
    chisq_path = open("C:\Users\Heathro\mg-gap\mg-gap\mg-gap-py\mg-gap\support_files\chisq.txt", "rU")
except:
    print("Error opening file paths")


# STEP 2 ----------
#  - run the vcf parser to determine b value for a SNP window of 1 (s = 1) 
#  - creates a snp list file called B1_new.txt, which genwin will read later
start_time = time.time()
print("Starting B processing at :", start_time)

# TODO question: what is the 'N' in the arguments for VCF_Analyzer?...only takes 3
# arguments, but with the 'N' it has 4 arguments
write_results = open("B1_new.txt", "w+")
write_results.write("CHR\tBP\tB\n")
snp_list = VCF_Analyzer.SNP_list(1, vcf_path, chisq_path)
for snp in snp_list:
    write_results.write("" + snp.Chromosome + '\t' + snp.Basepair + '\t' + snp.B_standard + "\n")
elapsed_time = time.time() - start_time
print("B processing time: ", elapsed_time)


# STEP 3 ----------
#  - calls the R script, GenWin 
#  - this reads the B1_new.txt file that was created above
#  - GenWin creates a file called "splinewindows.txt"
start_time = time.time()
print("Beginning R execution at", start_time)
# NOTE for testing: enter correct path
# TODO need to test this R call. Make sure Rscript is running from correct directory
rscript_path = "C:\Users\Heathro\mg-gap\mg-gap\mg-gap-py\mg-gap\support_files\GenWin_script_12_29_2016.R"
subprocess.call(["Rscript", rscript_path])
elapsed_time = time.time() - start_time
print("R exited successfully.\nRun time: ", elapsed_time)
 

# STEP 4 ----------
#  - calculate the median from the "splinewindows.txt" file created by GenWin
#  - run the b processing at that window and then b* processing 
#  - then move to java program?
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
splinepath = "N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/bin/Debug/splinewindows.txt";
# TODO add some error checking code, like try/catch?
print("GenWin file found, obtaining median window value...")
windows = []
with open(splinepath, "r") as sp:
    sp.readline() # read past the header line
    for line in sp:
        if ("CHRcol" not in line):
            cols = []
            cols.append(line.replace("\n", "").split("\t"))
            windows.append(int(cols[3]))
if len(windows) == 0:
    print("No window sizes to calculate!")
else:
     windows.sort()
     middle = len(windows) // 2
     median = float(windows[middle]) if (len(windows) != 0) else (float(windows[middle]) + float(windows[middle - 1])) / 2.0
     # TODO question: why are we making median a double/float if we are casting it
     # back into an integer in the next step


# STEP 5 ----------
#  - feed the median value back through b processing, then b* 
#  - need to hold on to the data for FDR 
#  - make a tab .csv
if median > 0:
    snpList = VCF_Analyzer.SNP_list(int(median), vcf_path, chisq_path)
    print("Re-analyzing for B* based on median window size ", median, " @ ", time.ctime)

    start_time = time.time()
    print("Starting FDR sorting processing at :", start_time)

    # Start the FDR process
    # This should be a user config variable in the future if they would prefer any different settings
    # TODO maybe create user input for fdr_input value
    fdr_input = 0.05
    print("Running FDR analysis at ", fdr_input, "...")
    
    fdrlist = FDR.Process(snpList, fdr_input)

    # Save this fdrlist to a "tab" separated values, use a .csv, file with headers 
    # TODO this is the same code as in step 2 - good candidate for a write snp list method
    write_results = open("FDR_analysis.txt", "w+")
    write_results.write("CHR\tBP\tB\n")
    for snp in fdrlist:
        write_results.write("" + snp.Chromosome + '\t' + snp.Basepair + '\t' + snp.B_standard + "\n")
    elapsed_time = time.time() - start_time
    print("Sort FDR processing time: ", elapsed_time)

  
    # Start getting the annotations !! Can stop here for the purposes of 
    # this program. 

    # TODO Something that may be of interest later. Connect to pythosome site to get Gene, description, and RNA SEQ data.