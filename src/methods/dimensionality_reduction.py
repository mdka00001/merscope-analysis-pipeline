from src.methods.dependencies import *
from src.methods.base import BaseMethod
from src.output.outputDTO import OutputDTO


class Dimensional(BaseMethod):
    def run(self):
        """
        Perform dimensionality reduction using PCA and UMAP.
        """
        # Load the input data
        print(f"Loading data from {self.adata}")
        adata = sc.read_h5ad(self.adata)

        # dividing by volume instead
        print("normalize total")
        sc.pp.normalize_total(adata)
        print("log transform")
        sc.pp.log1p(adata)
        print("scale")
        sc.pp.scale(adata, max_value=10)


        print("PCA")
        sc.tl.pca(adata, svd_solver="arpack")
        print("neighbors")
        sc.pp.neighbors(adata, 
                        n_neighbors=self.n_neighbors, 
                        n_pcs=self.n_pcs)
        print("UMAP")
        sc.tl.umap(adata)

        print("tSNE")
        sc.tl.tsne(adata, n_pcs=self.n_pcs)

        output = OutputDTO()
        output.add_data(adata)

        #plotting
        plots = ["pca_variance", "umap", "pca", "tsne"]

        for plot in plots:
            if plot == "pca_variance":
                sc.pl.pca_variance_ratio(adata, log=True, show=False)
                fig = plt.gcf()
                output.add_plot("pca_variance", fig)
            elif plot == "umap":
                sc.pl.umap(adata, color=["n_genes_by_counts", "total_counts"], show=False)
                fig = plt.gcf()
                output.add_plot("umap", fig)
            elif plot == "pca":
                sc.pl.pca(adata, color=["n_genes_by_counts", "total_counts"], show=False)
                fig = plt.gcf()
                output.add_plot("pca", fig)
            elif plot == "tsne":
                sc.pl.tsne(adata, color=["n_genes_by_counts", "total_counts"], show=False)
                fig = plt.gcf()
                output.add_plot("tsne", fig)


        return output