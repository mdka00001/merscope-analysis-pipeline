from src.methods.dependencies import *
import hdf5plugin
class OutputDTO:
    def __init__(self, data=None, plots=None):
        """
        data: dict of key results (e.g., filtered counts, clusters, annotations)
        plots: dict of plot objects (e.g., matplotlib figures)
        """
        self.data = data if data is not None else {}
        self.plots = plots if plots is not None else {}

    def add_data(self, data):
        self.data = data

    def add_plot(self, value):
        self.plots = value

    def get_data(self):
        return self.data

    def get_plots(self):
        return self.plots
    
    def save_plot(self, filename: str = "plot", directory: str = "plots", file_format: str = "png", dpi: int = 300):
        """
        Saves a matplotlib figure to a specified directory with the given filename.

        Parameters:
            figure (matplotlib.figure.Figure): The plot object to save.
            filename (str): Desired filename without extension.
            directory (str): Directory to save the plot (default is 'plots').
            file_format (str): Image format (e.g., 'png', 'jpg', 'pdf').
            dpi (int): Resolution of the saved image.
        """
        if not os.path.exists(directory):
            os.makedirs(directory)

        path = os.path.join(directory, f"{filename}.{file_format}")
        self.plots.savefig(fr"{path}")
        print(f"Plot saved to: {path}")

    def save_data(self, filename: str = "adata", directory: str = "data"):
        self.get_data().write(fr"{filename}.h5ad")
        print(f"Data saved to: {filename}.h5ad")