from src.methods.dependencies import *
from src.methods.base import BaseMethod
from src.output.outputDTO import OutputDTO

class CreateScanpyObject(BaseMethod):
    def run(self):
        """
        Create a Scanpy object from the input cell by gene matrix and cell metadata.
        """
        print(f"Loading data from {self.input_cell_by_gene} and {self.input_cell_metadata}")
        cell_by_gene = pd.read_csv(self.input_cell_by_gene, index_col=0)
        meta_cell = pd.read_csv(self.input_cell_metadata, index_col=0)

        # Create a Scanpy object
        meta_cell['barcodeCount'] = cell_by_gene.sum(axis=1)
        meta_gene = pd.DataFrame(index=cell_by_gene.columns.tolist())

        # Drop blanks for single cell analysis
        keep_genes = [x for x in cell_by_gene.columns.tolist() if 'Blank' not in x]
        cell_by_gene = cell_by_gene[keep_genes]
        meta_gene = meta_gene.loc[keep_genes]
        meta_gene['expression'] = cell_by_gene.sum(axis=0)

        vizgen_genes = meta_gene.index.tolist()

        #load metagene
        
        df_ref_panel_ini = pd.read_excel(self.ref_marker_panel, index_col=0)
        df_ref_panel = df_ref_panel_ini.iloc[1:,:1]
        df_ref_panel.index.name = None
        df_ref_panel.columns = ['Function']

        
        marker_genes = df_ref_panel[df_ref_panel['Function'].str.contains('marker')].index.tolist()
        common_marker_genes = list(set(meta_gene.index.tolist()).intersection(marker_genes))
        meta_gene.loc[common_marker_genes, 'Markers'] = df_ref_panel.loc[common_marker_genes, 'Function']
        meta_gene['Markers'] = meta_gene['Markers'].apply(lambda x: 'N.A.' if 'marker' not in str(x) else x)
        meta_gene['Markers'].value_counts()


        min_count = 50
        keep_cells = meta_cell[meta_cell['barcodeCount'] >= min_count].index.tolist()
        meta_cell = meta_cell.loc[keep_cells]
        cell_by_gene = cell_by_gene.loc[keep_cells]

        # reset index of cell_by_gene and meta_cell
        new_index = range(len(meta_cell.index.tolist()))
        meta_cell.index = new_index
        cell_by_gene.index = new_index
        cell_by_gene.shape

        adata = ad.AnnData(X=cell_by_gene.values, obs=meta_cell, var=meta_gene)
        adata.obsm["spatial"] = meta_cell[["center_x", "center_y"]].values
        adata.var_names_make_unique()

        # Annotate mitochondrial genes
        print("Identifying mitochondrial genes...")
        adata.var["mt"] = adata.var_names.str.startswith("mt-")

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=(50, 100, 200, 300), inplace=True
        )

        output = OutputDTO()
        output.add_data(adata)

        print(adata)

        # Generate and store QC plots
        qc_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
        for metric in qc_metrics:
            sc.pl.violin(
                adata,
                keys=metric,
                jitter=0.4,
                show=False
            )
            fig = plt.gcf()
            output.add_plot(f"{metric}_violin", fig)

        return output
