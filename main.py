from src.input.input_args import *
from src.methods.dependencies import *
from src.methods.base import *
from src.methods.create_scanpy_object import *
from src.output.outputDTO import *


def main():
    # Parse input arguments
    args = parse_args()

    if args.command == "create_scanpy_object":

        # Create a Scanpy object
        create_scanpy_object = CreateScanpyObject(args.input_cell_by_gene, args.input_cell_metadata)
        output = create_scanpy_object.run()

        output.save_plot(filename="qc_metrics", directory="plots", file_format="png", dpi=300)
        print("Scanpy object created and saved successfully.")
    else:
        print("Invalid command. Please use 'create_scanpy_object'.")
if __name__ == "__main__":
    main()