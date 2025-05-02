from src.methods.dependencies import *

class BaseMethod:
    def __init__(self, input_cell_by_gene, input_cell_metadata):
        self.input_cell_by_gene = input_cell_by_gene
        self.input_cell_metadata = input_cell_metadata

        
