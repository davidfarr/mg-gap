"""
Class SNP

Description: Sets up fields for each instantiation of a SNP. Getters 
and setters are accessed by the instance of the class. 
No private variables.
"""
class SNP:

    # constructor, sets all initial instance field values
    def __init__(self):
        self.basepair = 0
        self.chromosome = 0
        self.b_standard = 0.00
        self.b_star = 0.00
        self.raw_p = 0.00
        self.adjusted_p = ""
        self.description = ""
        self.old_identifier = "" # snnfold_X
        self.c_variance = 0.00
        self.t_variance = 0.00
        self.transformed_c_variance = 0.0
        self.transformed_t_variance = 0.0
        self.threshold_value = 0.0
        self.gene = ""
        self.originalindex = 0
        self.rnaseqp = ""
        self.fdr_rank = 0
