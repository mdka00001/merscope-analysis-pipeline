import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Parse input arguments for the Merscope analysis pipeline.")
    subparsers = parser.add_subparsers(dest='command', required=True, help='Sub-command help')
    
    # Add arguments

    # create scanpy object
    parser_a = subparsers.add_parser('create_scanpy_object', help='Create a scanpy object from input files.')
    parser_a.add_argument('--input_cell_by_gene', type=str, required=True, help='csv file of cell by genes.')
    parser_a.add_argument('--input_cell_metadata', type=str, required=True, help='csv file of cell metadata.')


    #filter cells
    parser_b = subparsers.add_parser('filter_cells', help='Filter cells based on metadata.')
    parser_b.add_argument('--adata', type=str, required=True, help='annotation data (scanpy object).')
    parser_b.add_argument('--min_counts', type=int, required=True, help='Minimum counts for filtering counts.', default=50)
    parser_b.add_argument('--min_cells', type=int, required=True, help='Minimum genes for filtering cells.', default=10)
    parser_b.add_argument('--n_genes_by_counts', type=int, required=True, help='Maximum n_genes_by_counts threshold.')
    parser_b.add_argument('--total_counts', type=int, required=True, help='Maximum total_counts threshold.')



    #clustering
    parser_d = subparsers.add_parser('clustering', help='Perform clustering on the data.')
    parser_d.add_argument('--adata', type=str, required=True, help='annotation data (scanpy object).')
    parser_d.add_argument('--resolution', type=float, required=True, help='Resolution parameter for clustering.')
    parser_d.add_argument('--n_pcs', type=int, required=True, help='Number of principal components to compute.')
    parser_d.add_argument('--n_neighbors', type=int, required=True, help='Number of neighbors for UMAP.')

    #cluster annotation
    parser_e = subparsers.add_parser('cluster_annotation', help='Annotate clusters.')
    parser_e.add_argument('--adata', type=str, required=True, help='annotation data (scanpy object).')
    parser_e.add_argument('--ref_marker_panel', type=str, required=True, help='Reference marker genes dataset for cluster annotation.')

    #visualize spatial map
    parser_f = subparsers.add_parser('visualize_spatial_map', help='Visualize spatial map of the data.')
    parser_f.add_argument('--adata', type=str, required=True, help='annotation data (scanpy object).')
    #
    
    # Parse arguments
    args = parser.parse_args()
    
    return args