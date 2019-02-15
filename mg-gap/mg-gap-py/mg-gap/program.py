# import python packages
import time
import subprocess # for running an R script

# import classes
import VCF_Analyzer
import FDR 


# STEP 1 ----------
#  - set up the vcf and chisq filepaths
#  - NOTE for testing, enter the correct paths
vcf_path = "N:/app dev/scoville research/program files/dev migration for windows/vcf files/Ali_w_767.vcf";
chisq_path = "C:/Users/David/Documents/GitHub/mg-gap/mg-gap/mg-gap/support files/chisq.txt";


# STEP 2 ----------
#  - run the vcf parser to determine b value for a SNP window of 1 (s = 1) 
#  - creates a snp list file called B1_new.txt, which genwin will read later
start_time = time.time()
print("Starting B processing at :", start_time)

# TODO question: what is the 'N' in the arguments for VCF_Analyzer?...only takes 3
# arguments, but with the 'N' is 4
write_results = open("B1_new.txt", "w+")
write_results.write("CHR\tBP\tB\n")
for snp in VCF_Analyzer.SNP_list(1, vcf_path, 'N', chisq_path):
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
rscript_path = "N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/support files/GenWin_script_12_29_2016.R"
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
# TODO add some error checking code, like try/catch
print("GenWin file found, obtaining median window value...")
windows = []
with open(splinepath, "r") as sp:
    sp.readline() # read past the header
    for line in sp:
        if ("CHRcol" not in line):
            cols = []
            cols.append(line.replace("\n", "").split("\t"))
            windows.append()
            # NOTE !!! Stopped here


# STEP 5 ----------
#  - feed the median value back through b processing, then b* 
#  - need to hold on to the data for FDR though
median = 6
if median > 0:
    # TODO need to find equivalent "Convert" method in python see commented code below
    snpList = VCF_Analyzer.SNP_list(Convert.ToInt32(median), vcf_path, chisq_path)
    print("Re-analyzing for B* based on median window size ", median, " @ ", time.ctime)

    # Start the FDR process
    # This should be a user config variable in the future if they would prefer any different settings
    # TODO maybe create user input for fdr_input value
    fdr_input = 0.05
    print("Running FDR analysis at ", fdr_input, "...")
    # TODO save this fdrlist to a "tab" separated values, use a .csv, file with headers 
    fdrlist = FDR.Process(snpList, fdr_input)
    
    # Start getting the annotations !! May stop here for the purposes of 
    # this program. TODO ask about this. Yes we canstop here.

    # TODO Something that may be of interest later. Connect to pythosome site to get Gene, description, and RNA SEQ data.


# TODO Look up "Convert" method equivalent in python
"""
namespace System
{
    //
    // Summary:
    //     Converts a base data type to another base data type.
    public static class Convert {

        public static int ToInt32(int value);
        //
        // Summary:
        //     Converts the value of the specified 64-bit signed integer to an equivalent 32-bit
        //     signed integer.
        //
        // Parameters:
        //   value:
        //     The 64-bit signed integer to convert.
        //
        // Returns:
        //     A 32-bit signed integer that is equivalent to value.
        //
        // Exceptions:
        //   T:System.OverflowException:
        //     value is greater than System.Int32.MaxValue or less than System.Int32.MinValue.
}
"""
