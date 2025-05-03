from src.methods.dependencies import *
import hdf5plugin

class OutputDTO:
    def __init__(self, data=None, plots=None):
        """
        data: AnnData or similar primary output
        plots: dict of named matplotlib figures (e.g., {'qc': fig1, 'umap': fig2})
        """
        self.data = data
        self.plots = plots if plots is not None else {}

    def add_data(self, data):
        self.data = data

    def add_plot(self, name: str, figure):
        """Add a named plot to the collection."""
        self.plots[name] = figure

    def get_data(self):
        return self.data

    def get_plots(self):
        return self.plots

    def save_plot(self, directory: str = "plots", file_format: str = "png", dpi: int = 300):
        """
        Saves all plots in the `plots` dictionary to the specified directory.
        """
        if not os.path.exists(directory):
            os.makedirs(directory)

        for name, fig in self.plots.items():
            path = os.path.join(directory, f"{name}.{file_format}")
            fig.savefig(path, dpi=dpi)
            print(f"Plot saved to: {path}")

    def save_data(self, filename: str = "adata", directory: str = "data"):
        """Saves the AnnData object to an h5ad file."""
        if not os.path.exists(directory):
            os.makedirs(directory)

        full_path = os.path.join(directory, f"{filename}.h5ad")
        self.data.write(full_path, compression="gzip")
        print(f"Data saved to: {full_path}")
