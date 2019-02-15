# import python packages
import time

# import classes
import VCF_Analyzer
import FDR 

# set up the vcf and chisq filepaths
# TODO enter the correct paths
vcf_path = "N:/app dev/scoville research/program files/dev migration for windows/vcf files/Ali_w_767.vcf";
chisq_path = "C:/Users/David/Documents/GitHub/mg-gap/mg-gap/mg-gap/support files/chisq.txt";

# run the vcf parser for SNP window of 1
start_time = time.time()
print("Starting B processing at :", start_time)

# TODO Maybe pass variable instead of writing to file, the variable is a list of SNPs
# then again, check what genwin needs, it may NEED a file to read in the script

# open a new file, w is write, + is create if not already in directory
# TODO question: what is the 'N' in the arguments for VCF_Analyzer?...only takes 3
# arguments, but with the 'N' is 4
write_results = open("B1_new.txt", "w+")
write_results.write("CHR\tBP\tB\n")
for snp in VCF_Analyzer.SNP_list(1, vcf_path, 'N', chisq_path):
    write_results.write("" + snp.Chromosome + '\t' + snp.Basepair + '\t' + snp.B_standard + "\n")
elapsed_time = time.time() - start_time
print("B processing time: ", elapsed_time)

# after b values in a window of s = 1 is determined above,
# this next block feeds results into genwin
start_time = time.time()
print("Beginning R execution at", start_time)
# create R engine, do try/ catch?
print("Successfully created R engine instance. Evaluating script...")
#


# feed the median value back through b processing, then b* 
# need to hold on to the data for FDR though
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
