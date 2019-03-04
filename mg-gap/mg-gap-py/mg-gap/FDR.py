"""
Class FDR

Description: 
The FDR process is:
    1. Sort B* list by p-value
    2. Add new property FDR where FDR = 0.1 * index of the ranked SNP / (count of SNPs in the file / 2)
    3*. FDR can be changed... 0.1 above is FDR of 10 and 0.05 is FDR 5

Editing thoughts for later: this class contains only one method.
No fields exist for this class when it is instantiated. Perhaps 
move method to the main program.py file as a method outside of main?
Need to look up python equivalents of C# list operations: 
    .OrderBy()
    .ToList()
    .Count()
"""
# class imports
# may need to import snp class?

#class FDR:
def process(bs_list, fdr_selected):
       
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
    print("\nSignificant B* (the min B* after removing SNPs over threshold) %s" % min(sorted_list, key = lambda snp: snp.b_star))

    return sorted_list # return a list of SNPs