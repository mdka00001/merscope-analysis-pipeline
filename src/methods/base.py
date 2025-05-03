from src.methods.dependencies import *

class BaseMethod:
    def __init__(self, input_cell_by_gene=None, input_cell_metadata=None,
                 adata = None, min_counts=50, min_cells=10,
                 n_genes_by_counts=None, total_counts=None,
                 n_pcs=None, n_neighbors=None,
                 resolution=None, ref_marker_panel=None, tsne=None):
        
        self.input_cell_by_gene = input_cell_by_gene if input_cell_by_gene is not None else {}
        self.input_cell_metadata = input_cell_metadata if input_cell_metadata is not None else {}
        self.adata = adata if adata is not None else {}
        self.min_counts = min_counts
        self.min_cells = min_cells
        self.n_genes_by_counts = n_genes_by_counts if n_genes_by_counts is not None else {}
        self.total_counts = total_counts if total_counts is not None else {}
        self.n_pcs = n_pcs if n_pcs is not None else {}
        self.n_neighbors = n_neighbors if n_neighbors is not None else {}
        self.resolution = resolution if resolution is not None else {}
        self.ref_marker_panel = ref_marker_panel if ref_marker_panel is not None else {}
        self.tsne = tsne if tsne is not None else {}


        
