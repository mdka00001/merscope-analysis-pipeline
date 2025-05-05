from src.input.input_args import *
from src.methods.dependencies import *
from src.methods.base import *
from src.methods.create_scanpy_object import *
from src.output.outputDTO import *
from src.methods.filter_cells import *
from src.methods.clustering import *
from src.methods.cluster_annotation import *
from src.methods.visualizer import *

def main():
    # Parse input arguments
    args = parse_args()

    if args.command == "create_scanpy_object":

        print("Creating Scanpy object...")
        print("This may take a while depending on the size of the input data.")

        # Create a Scanpy object
        create_scanpy_object = CreateScanpyObject(
            input_cell_by_gene=args.input_cell_by_gene, 
            input_cell_metadata=args.input_cell_metadata,
            ref_marker_panel=args.ref_marker_panel)
        
        
        output = create_scanpy_object.run()

        output.save_plot(directory="plots", file_format="png", dpi=300)
        print("Scanpy object created and saved successfully.")

        output.save_data(filename="adata", directory="data")

    elif args.command == "filter_cells":
        print("Filtering cells...")
        filter_cells = FilterCells(
            adata=args.adata,
            min_counts=args.min_counts,
            min_cells=args.min_cells,
            n_genes_by_counts=args.n_genes_by_counts,
            total_counts=args.total_counts
        )

        output = filter_cells.run()

        output.save_plot(directory="plots", file_format="png", dpi=300)
        output.save_data(filename="filtered_adata", directory="data")
        print("Cells filtered successfully.")

    elif args.command == "clustering":
        print("Performing clustering...")
        clustering = Clustering(
            adata=args.adata,
            resolution=args.resolution,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            tsne=args.tsne
        )

        output = clustering.run()
        output.save_plot(directory="plots", file_format="png", dpi=300)
        print("Clustering completed successfully.")

        output.save_data(filename="clustered_adata", directory="data")

    elif args.command == "cluster_annotation":
        print("Annotating clusters...")
        cluster_annotation = ClusterAnnotation(
            adata=args.adata,
            input_cell_by_gene=args.input_cell_by_gene
        )

        output = cluster_annotation.run()
        output.save_plot(directory="plots", file_format="png", dpi=300)
        print("Cluster annotation completed successfully.")

        output.save_data(filename="annotated_adata", directory="data")
    elif args.command == "visualize_spatial_map":
        print("Visualizing spatial map...")
        visualize_spatial_map = Visualizer(
            adata=args.adata,
            ref_marker_panel=args.ref_marker_panel,
            input_cell_by_gene=args.input_cell_by_gene,
            input_cell_metadata=args.input_cell_metadata
        )

        visualize_spatial_map.run()


    
    else:
        print("Invalid command. Please use 'create_scanpy_object'.")
if __name__ == "__main__":
    main()