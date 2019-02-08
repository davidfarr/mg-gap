"""
Class SNP

Description: Sets up fields for each instantiation of a SNP. Getters 
and setters are accessed by the instance of the class. 
No private variables.
"""
class SNP:

    # constructor, sets all initial instance field values
    def __init__(self):
        basepair = 0
        chromosome = 0
        b_standard = 0.00
        b_star = 0.00
        raw_p = 0.00
        adjusted_p = ""
        description = ""
        old_identifier = "" # snnfold_X
        c_variance = 0.00
        t_variance = 0.00
        transformed_c_variance = 0.0
        transformed_t_variance = 0.0
        threshold_value = 0.0
        gene = ""
        originalindex = 0
        rnaseqp = ""
        fdr_rank = 0
