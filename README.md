# Merscope Analysis Pipeline

This repository provides a small command line pipeline for analysing Merscope datasets with [Scanpy](https://scanpy.readthedocs.io/) and related tools.  The code allows you to build a `scanpy` object from CSV inputs, filter cells and genes, and run dimensionality reduction with PCA, UMAP and t‑SNE.

## Directory Structure

```
.
├── main.py               # Entry point for the CLI
├── src
│   ├── input
│   │   └── input_args.py           # Argument parser and sub‑commands
│   ├── methods
│   │   ├── base.py                 # Shared options for pipeline steps
│   │   ├── create_scanpy_object.py # Build AnnData object and QC plots
│   │   ├── dimensionality_reduction.py
│   │   ├── filter_cells.py         # Basic cell/feature filtering
│   │   └── dependencies.py         # Common imports used by the methods
│   └── output
│       └── outputDTO.py            # Helper class for data/plot output
```

## Usage

The pipeline is executed via `main.py` and exposes several sub‑commands:

### 1. Create a Scanpy object

```
python main.py create_scanpy_object \
  --input_cell_by_gene <cell_by_gene.csv> \
  --input_cell_metadata <cell_metadata.csv>
```

This step reads the CSV files, drops blank genes, calculates QC metrics (including mitochondrial content) and produces violin plots for metrics such as `n_genes_by_counts` and `total_counts`. The result is stored as an AnnData object.

### 2. Filter cells

```
python main.py filter_cells \
  --adata <adata.h5ad> \
  --min_counts 50 \
  --min_cells 10 \
  --n_genes_by_counts 4000 \
  --total_counts 40000
```

Filters cells and genes using Scanpy utility functions and simple thresholds. The filtered AnnData object is saved to disk.

### 3. Dimensionality reduction

```
python main.py dimensionality_reduction \
  --adata <filtered_adata.h5ad> \
  --n_pcs 50 \
  --n_neighbors 15
```

Normalises the data, logs, scales and then computes PCA, UMAP and t‑SNE embeddings. Plots for each step are written to the `plots` directory.

## OutputDTO

Results from each method are wrapped in an `OutputDTO` object located in `src/output/outputDTO.py`. It stores the AnnData result and any generated matplotlib figures, and provides utility methods `save_plot()` and `save_data()` to persist them.

## Dependencies

Dependencies are imported in `src/methods/dependencies.py` and include common scientific libraries (`numpy`, `pandas`, `scanpy`, `squidpy`, etc.). Some optional libraries like `clustergrammer2` are also imported.

## Notes

Only the three sub‑commands shown above are implemented in `main.py`. Additional commands defined in `input_args.py` (`clustering`, `cluster_annotation`, `visualize_spatial_map`) are present but not yet implemented.

