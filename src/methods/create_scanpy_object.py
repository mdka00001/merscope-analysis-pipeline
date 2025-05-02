from src.methods.dependencies import *
from src.methods.base import BaseMethod
from src.output.outputDTO import OutputDTO

class CreateScanpyObject(BaseMethod):
    def run(self):
        """
        Create a Scanpy object from the input cell by gene matrix and cell metadata.
        """

        # Load the input data
        print(f"Loading data from {self.input_cell_by_gene} and {self.input_cell_metadata}")
        cell_by_gene = pd.read_csv(self.input_cell_by_gene, index_col=0)
        meta_cell = pd.read_csv(self.input_cell_metadata, index_col=0)

        # Create a Scanpy object
        meta_cell['barcodeCount'] = cell_by_gene.sum(axis=1)

        # initialize meta_gene
        meta_gene = pd.DataFrame(index=cell_by_gene.columns.tolist())

        # drop blanks for single cell analysis
        keep_genes = [x for x in cell_by_gene.columns.tolist() if 'Blank' not in x]

        cell_by_gene = cell_by_gene[keep_genes]
        meta_gene = meta_gene.loc[keep_genes]

        meta_gene['expression'] = cell_by_gene.sum(axis=0)

        vizgen_genes = meta_gene.index.tolist()

        adata = sc.AnnData(X=cell_by_gene.values, obs=meta_cell, var=meta_gene)
    

        #unique gene name
        adata.var_names_make_unique()

        ##mitochondrial genes
        print(f"mitochondrial genes")
        adata.var["mt"] = adata.var_names.str.startswith("mt-")


        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=(50, 100, 200, 300), inplace=True
        )


        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False
        )

        fig = plt.gcf()

        output = OutputDTO()
        output.add_plot(fig)
        output.add_data(adata)

        return output
