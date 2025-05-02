from src.methods.dependencies import *
from src.methods.base import BaseMethod
from src.output.outputDTO import OutputDTO


class FilterCells(BaseMethod):
    def run(self):
        """
        Filter cells based on the number of genes expressed and total counts.
        """
        # Load the input data
        print(f"Loading data from {self.adata}")
        adata = sc.read_h5ad(self.adata)
        
        sc.pp.filter_cells(adata, min_counts=self.min_counts)
        sc.pp.filter_genes(adata, min_cells=self.min_cells)

        # Filter cells based on n_genes_by_counts and total_counts
        adata = adata[(adata.obs["n_genes_by_counts"] < self.n_genes_by_counts) ].copy()
        # Filter cells based on n_genes_by_counts and total_counts
        adata = adata[(adata.obs["total_counts"] < self.total_counts) ].copy()

        output = OutputDTO()
        output.add_data(adata)
        return output