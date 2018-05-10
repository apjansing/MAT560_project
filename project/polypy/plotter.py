import pylab as pl


class plotter:
    """docstring for plotter"""

    def __init__(self,
                 figure=pl.figure(),
                 rows=1,
                 columns=1,
                 xlims=[-10, 10],
                 ylims=[-10, 10],
                 tight_layout_pad=2):
        self.fig = figure
        self.rows = rows
        self.columns = columns
        self.xlims = xlims
        self.ylims = ylims
        pl.tight_layout(tight_layout_pad)

    def createSubplot(self, to_plot, loc, title, x_lims=None, y_lims=None):
        x_lims = self.xlims if x_lims == None else x_lims
        y_lims = self.ylims if y_lims == None else y_lims
        g = pl.subplot(self.rows, self.columns, loc)
        g.set_title(title)
        for points in to_plot:
            g.plot(points[0], points[1])
        g.set_xlim(x_lims)
        g.set_ylim(y_lims)

    def setLayout(self, pad=2):
        pl.tight_layout(pad)

    def show(self):
        pl.show()