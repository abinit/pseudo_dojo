import numpy as np
import matplotlib.pyplot as plt

from matplotlib.collections import PatchCollection
from ptplotter.plotter import ElementDataPlotter


class ElementDataPlotterRangefixer(ElementDataPlotter):
    """
    modify plotter to allow setting the clim for the plot
    """

    def draw(self, colorbars=True, **kwargs):
        """Replace super().draw."""
        self.cbars = []
        clims = kwargs.get('clims', None)
        n = len(self.collections)
        if clims is None:
            clims = [None]*n
        elif len(clims) == 1:
            clims = [clims[0]]*n
        elif len(clims) == n:
            pass
        else:
            raise RuntimeError('incorrect number of clims provided in draw')

        for coll, cmap, label, clim in zip(self.collections, self.cmaps, self.cbar_labels, clims):
            #print(clim)
            pc = PatchCollection(coll, cmap=cmap)
            pc.set_clim(vmin=clim[0],vmax=clim[1])
            #print(pc.get_clim())
            pc.set_array(np.array([p.value for p in coll]))
            self._ax.add_collection(pc)

            if colorbars:
                options = {'orientation':'horizontal', 'pad':0.05, 'aspect':60,}
                options.update(kwargs.get('colorbar-options', {}))
                cbar = plt.colorbar(pc, **options)
                cbar.set_label(label)
                self.cbars.append(cbar)
        fontdict = kwargs.get('font', {'color':'white'})

        for s in self.squares:
            if not s.label: continue
            x = s.x + s.dx/2
            y = s.y + s.dy/2
            self._ax.text(x, y, s.label, ha='center', va='center', fontdict=fontdict)

        qs_labels = [k.split('[')[0] for k in self.labels]

        if self.guide_square:
            self.guide_square.set_labels(qs_labels)
            pc = PatchCollection(self.guide_square.patches, match_original=True)
            self._ax.add_collection(pc)
        self._ax.autoscale_view()
