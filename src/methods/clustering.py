from src.methods.dependencies import *
from src.methods.base import BaseMethod
from src.output.outputDTO import OutputDTO


class Clustering(BaseMethod):
    def run(self):
        """
        Perform clustering on the input data using Louvain and Leiden algorithms.
        """
        # Load the input data
        print(f"Loading data from {self.adata}")
        adata = sc.read_h5ad(self.adata)

        print("neighbors")
        sc.pp.neighbors(adata, 
                        n_neighbors=self.n_neighbors, 
                        n_pcs=self.n_pcs)
        print("UMAP")
        sc.tl.umap(adata)

        print("tSNE")
        sc.tl.tsne(adata, n_pcs=self.n_pcs)

        print("Leiden")
        sc.tl.leiden(adata, resolution=self.resolution)

        output = OutputDTO()
        output.add_data(adata)

        # Plotting
        plots = ["umap_leiden", "tSNE_leiden"]

        for plot in plots:
            if plot == "tSNE_leiden":
                sc.pl.tsne(adata, color=["leiden"], show=False)
                fig = plt.gcf()
                output.add_plot("tSNE_leiden", fig)
            elif plot == "umap_leiden":
                sc.pl.umap(adata, color=["leiden"], show=False)
                fig = plt.gcf()
                output.add_plot("umap_leiden", fig)

        return output