from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import abc

from library import file_utils
from library.graphs.graph_utils import legend_patches_and_labels


class GraphBase(metaclass=abc.ABCMeta):

    def __init__(self, **kwargs):
        self.title = kwargs.get("title")
        self.x_label = kwargs.get("x_label")
        self.y_label = kwargs.get("y_label")
        self.legend = kwargs.get('legend')
        self.plot_methods = [self.plot]
        self.figsize = None

    def decorations(self, ax):
        if self.title:
            ax.set_title(self.title)
        if self.x_label:
            ax.set_xlabel(self.x_label)
        if self.y_label:
            ax.set_ylabel(self.y_label)

    @abc.abstractmethod
    def plot(self, ax):
        pass

    def post_plot(self, ax):
        if self.legend:
            patches, labels = legend_patches_and_labels(self.legend)

            ax.legend(patches, labels, loc='upper left')

    def figure(self, figure):
        """ a hook method if you want to do something about the figure """
        pass

    def plot_figure(self, figure):
        ax = figure.add_subplot(1, 1, 1)

        self.decorations(ax)
        for plot in self.plot_methods:
            plot(ax)  # Implementation
        self.post_plot(ax)

    def save(self, filename_or_obj, figsize=None):
        figsize = figsize or self.figsize

        figure = Figure(figsize=figsize)
        figure.patch.set_facecolor('white')

        self.plot_figure(figure)
        self.figure(figure)

        if isinstance(filename_or_obj, str):
            file_utils.mk_path_for_file(filename_or_obj)

        canvas = FigureCanvasAgg(figure)
        canvas.print_png(filename_or_obj)
