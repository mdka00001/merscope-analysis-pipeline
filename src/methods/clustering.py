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

        print("Leiden")
        sc.tl.leiden(adata, resolution=self.resolution)

        output = OutputDTO()
        output.add_data(adata)

        # Plotting
        plots = ["umap_leiden"]

        for plot in plots:
            if plot == "umap_louvain":
                sc.pl.umap(adata, color=["louvain"], show=False)
                fig = plt.gcf()
                output.add_plot("umap_louvain", fig)
            elif plot == "umap_leiden":
                sc.pl.umap(adata, color=["leiden"], show=False)
                fig = plt.gcf()
                output.add_plot("umap_leiden", fig)

        return output