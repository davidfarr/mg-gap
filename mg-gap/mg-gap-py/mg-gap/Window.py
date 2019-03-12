"""
Class Window

Description: Keeps track of windows to find median window value to analyze DNA

"""

class Window:
    # constructor, sets all initial instance field values
    def __init__(self):
        self.chromosome = 0
        self.snps = [] # list of SNPs
        self.window_id = 0
        self.window_size = 0

