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



        # dividing by volume instead
        print("normalize total")
        sc.pp.normalize_total(adata)
        print("log transform")
        sc.pp.log1p(adata)
        print("scale")
        sc.pp.scale(adata, max_value=10)


        print("PCA")
        sc.tl.pca(adata, svd_solver="arpack")


        output = OutputDTO()
        output.add_data(adata)

        plots = ["pca_variance", "pca"]

        for plot in plots:
            if plot == "pca_variance":
                sc.pl.pca_variance_ratio(adata, log=True, show=False)
                fig = plt.gcf()
                output.add_plot("pca_variance", fig)
            elif plot == "pca":
                sc.pl.pca(adata, color=["n_genes_by_counts", "total_counts"], show=False)
                fig = plt.gcf()
                output.add_plot("pca", fig)
            

        return output