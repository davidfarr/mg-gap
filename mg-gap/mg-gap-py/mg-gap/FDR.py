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
import SNP

class FDR:
    def process(bs_list, fdr_selected):
        # TODO add all the code :)
       
        sortedList = [] # place holder variable
        return sortedList # return a list of SNPs