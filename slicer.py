"""Hacked version of @prisae's 3D slicer in Discretize to interactive slice a
3D PyVista mesh using Matplotlib"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
import matplotlib.cm as cmx

import pyvista as pv
import numpy as np





class Slicer(object):
    """Plot slices of a 3D volume, interactively (scroll wheel).
    If called from a notebook, make sure to set
        %matplotlib notebook
    The straight forward usage for the Slicer is through, e.g., a
    `TensorMesh`-mesh, by accessing its `mesh.plot_3d_slicer`.
    If you, however, call this class directly, you have first to initiate a
    figure, and afterwards connect it:
    >>> # You have to initialize a figure
    >>> fig = plt.figure()
    >>> # Then you have to get the tracker from the Slicer
    >>> tracker = discretize.View.Slicer(mesh, Lpout)
    >>> # Finally you have to connect the tracker to the figure
    >>> fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
    >>> plt.show()
    **Parameters**
    v : array
        Data array of length self.nC.
    xslice, yslice, zslice : floats, optional
        Initial slice locations (in meter);
        defaults to the middle of the volume.
    vType: str
        Type of visualization. Default is 'CC'.
        One of ['CC', 'Fx', 'Fy', 'Fz', 'Ex', 'Ey', 'Ez'].
    view : str
        Which component to show. Defaults to 'real'.
        One of  ['real', 'imag', 'abs'].
    axis : 'xy' (default) or 'yx'
        'xy': horizontal axis is x, vertical axis is y. Reversed otherwise.
    transparent : 'slider' or list of floats or pairs of floats, optional
        Values to be removed. E.g. air, water.
        If single value, only exact matches are removed. Pairs are treated as
        ranges. E.g. [0.3, [1, 4], [-np.infty, -10]] removes all values equal
        to 0.3, all values between 1 and 4, and all values smaller than -10.
        If 'slider' is provided it will plot an interactive slider to choose
        the shown range.
    clim : None or list of [min, max]
        For pcolormesh (vmin, vmax).
    xlim, ylim, zlim : None or list of [min, max]
        Axis limits.
    aspect : 'auto', 'equal', or num
        Aspect ratio of subplots. Defaults to 'auto'.
        A list of two values can be provided. The first will be for the
        XY-plot, the second for the XZ- and YZ-plots, e.g. ['equal', 2] to have
        the vertical dimension exaggerated by a factor of 2.
        WARNING: For anything else than 'auto', unexpected things might happen
                 when zooming, and the subplot-arrangement won't look pretty.
    grid : list of 3 int
        Number of cells occupied by x, y, and z dimension on plt.subplot2grid.
    pcolorOpts : dictionary
        Passed to pcolormesh.
    """

    def __init__(self, mesh, scalars, xslice=None, yslice=None, zslice=None,
                 axis='xy', transparent=None,
                 clim=None, xlim=None, ylim=None, zlim=None, aspect='auto',
                 grid=[2, 2, 1], pcolorOpts=None):
        """Initialize interactive figure."""

        # 0. Some checks, not very extensive

        # (a) Mesh dimensionality
        # assumes 3D for PyVista.... not always be safe

        # (b) vType
        # Eh, I'm not sure how to plot cell data in MPL


        # (c) vOpts  # Not yet working for 'vec'


        # 1. Store relevant data

        # Store data in self
        self.mesh = mesh.ctp()

        # Axis: Default ('xy'): horizontal axis is x, vertical axis is y.
        # Reversed otherwise.
        self.yx = axis == 'yx'

        if xslice == None:
            xslice = mesh.center[0]
        if yslice == None:
            yslice = mesh.center[1]
        if zslice == None:
            zslice = mesh.center[2]

        self.xind = xslice
        self.yind = yslice
        self.zind = zslice

        # Aspect ratio
        if isinstance(aspect, (list, tuple)):
            aspect1 = aspect[0]
            aspect2 = aspect[1]
        else:
            aspect1 = aspect
            aspect2 = aspect
        if aspect2 in ['auto', 'equal']:
            aspect3 = aspect2
        else:
            aspect3 = 1.0/aspect2

        self.mesh.set_active_scalar(scalars)
        # Store min and max of all dat
        if clim is None:
            clim = self.mesh.get_data_range()
        self.pc_props = {'vmin': clim[0], 'vmax': clim[1]}

        # 2. Start populating figure

        # Get plot2grid dimension
        figgrid = (grid[0]+grid[2], grid[1]+grid[2])

        # Create subplots
        self.fig = plt.gcf()
        self.fig.subplots_adjust(wspace=.075, hspace=.1)

        # X-Y
        self.ax1 = plt.subplot2grid(figgrid, (0, 0), colspan=grid[1],
                                    rowspan=grid[0], aspect=aspect1)
        if self.yx:
            self.ax1.set_ylabel('x')
            if ylim is not None:
                self.ax1.set_xlim([ylim[0], ylim[1]])
            if xlim is not None:
                self.ax1.set_ylim([xlim[0], xlim[1]])
        else:
            self.ax1.set_ylabel('y')
            if xlim is not None:
                self.ax1.set_xlim([xlim[0], xlim[1]])
            if ylim is not None:
                self.ax1.set_ylim([ylim[0], ylim[1]])
        self.ax1.xaxis.set_ticks_position('top')
        plt.setp(self.ax1.get_xticklabels(), visible=False)

        # X-Z
        self.ax2 = plt.subplot2grid(figgrid, (grid[0], 0), colspan=grid[1],
                                    rowspan=grid[2], sharex=self.ax1,
                                    aspect=aspect2)
        self.ax2.yaxis.set_ticks_position('both')
        if self.yx:
            self.ax2.set_xlabel('y')
            if ylim is not None:
                self.ax2.set_xlim([ylim[0], ylim[1]])
        else:
            self.ax2.set_xlabel('x')
            if xlim is not None:
                self.ax2.set_xlim([xlim[0], xlim[1]])
        self.ax2.set_ylabel('z')
        if zlim is not None:
            self.ax2.set_ylim([zlim[0], zlim[1]])

        # Z-Y
        self.ax3 = plt.subplot2grid(figgrid, (0, grid[1]), colspan=grid[2],
                                    rowspan=grid[0], sharey=self.ax1,
                                    aspect=aspect3)
        self.ax3.yaxis.set_ticks_position('right')
        self.ax3.xaxis.set_ticks_position('both')
        self.ax3.invert_xaxis()
        plt.setp(self.ax3.get_yticklabels(), visible=False)
        if self.yx:
            if xlim is not None:
                self.ax3.set_ylim([xlim[0], xlim[1]])
        else:
            if ylim is not None:
                self.ax3.set_ylim([ylim[0], ylim[1]])
        if zlim is not None:
            self.ax3.set_xlim([zlim[1], zlim[0]])

        # Cross-line properties
        # We have to lines, a thick white one, and in the middle a thin black
        # one, to assure that the lines can be seen on dark and on bright
        # spots.
        self.clpropsw = {'c': 'w', 'lw': 2, 'zorder': 10}
        self.clpropsk = {'c': 'k', 'lw': 1, 'zorder': 11}

        # Add pcolorOpts
        if pcolorOpts is not None:
            self.pc_props.update(pcolorOpts)

        # Initial draw
        self.update_xy()
        self.update_xz()
        self.update_zy()

        # Create colorbar
        plt.colorbar(self.zy_pc, pad=0.15)

#         # Remove transparent value
#         if isinstance(transparent, str) and transparent.lower() == 'slider':
#             # Sliders
#             self.ax_smin = plt.axes([0.7, 0.11, 0.15, 0.03])
#             self.ax_smax = plt.axes([0.7, 0.15, 0.15, 0.03])

#             # Limits slightly below/above actual limits, clips otherwise
#             self.smin = Slider(self.ax_smin, 'Min', *clim, valinit=clim[0])
#             self.smax = Slider(self.ax_smax, 'Max', *clim, valinit=clim[1])

#             def update(val):
#                 self.v.mask = False  # Re-set
#                 self.v = np.ma.masked_outside(self.v.data, self.smin.val,
#                                               self.smax.val)
#                 # Update plots
#                 self.update_xy()
#                 self.update_xz()
#                 self.update_zy()

#             self.smax.on_changed(update)
#             self.smin.on_changed(update)

#         elif transparent is not None:

#             # Loop over values
#             for value in transparent:
#                 # If value is a list/tuple, we treat is as a range
#                 if isinstance(value, (list, tuple)):
#                     self.v = np.ma.masked_inside(self.v, value[0], value[1])
#                 else: # Exact value
#                     self.v = np.ma.masked_equal(self.v, value)

#             # Update plots
#             self.update_xy()
#             self.update_xz()
#             self.update_zy()

        # 3. Keep depth in X-Z and Z-Y in sync

        def do_adjust():
            """Return True if z-axis in X-Z and Z-Y are different."""
            one = np.array(self.ax2.get_ylim())
            two = np.array(self.ax3.get_xlim())[::-1]
            return sum(abs(one - two)) > 0.001  # Difference at least 1 m.

        def on_ylims_changed(ax):
            """Adjust Z-Y if X-Z changed."""
            if do_adjust():
                self.ax3.set_xlim([self.ax2.get_ylim()[1],
                                   self.ax2.get_ylim()[0]])

        def on_xlims_changed(ax):
            """Adjust X-Z if Z-Y changed."""
            if do_adjust():
                self.ax2.set_ylim([self.ax3.get_xlim()[1],
                                   self.ax3.get_xlim()[0]])

        self.ax3.callbacks.connect('xlim_changed', on_xlims_changed)
        self.ax2.callbacks.connect('ylim_changed', on_ylims_changed)

        self.fig.canvas.mpl_connect('scroll_event', self.onscroll)



    def onscroll(self, event):
        """Update index and data when scrolling."""
        print("onscroll")

        # Get scroll direction
        if event.button == 'up':
            pm = 1
        else:
            pm = -1

        # Update slice index depending on subplot over which mouse is
        if event.inaxes == self.ax1:    # X-Y
            self.zind = (self.zind + pm)
            self.update_xy()
        elif event.inaxes == self.ax2:  # X-Z
            if self.yx:
                self.xind = (self.xind + pm)
            else:
                self.yind = (self.yind + pm)
            self.update_xz()
        elif event.inaxes == self.ax3:  # Z-Y
            if self.yx:
                self.yind = (self.yind + pm)
            else:
                self.xind = (self.xind + pm)
            self.update_zy()

        plt.draw()



    def update_xy(self):
        """Update plot for change in Z-index."""
        print("updating xy")
        # Clean up
        self._clear_elements(['xy_pc', 'xz_ahw', 'xz_ahk', 'zy_avw', 'zy_avk'])

        # Draw X-Y slice
        cntr = self.mesh.center
        cntr[2] = self.zind
        self.xy_slc = self.mesh.slice(normal="z", origin=cntr, generate_triangles=True)
        pts = self.xy_slc.points
        tri = self.xy_slc.faces.reshape((-1,4))[:, 1:]
        val = self.xy_slc.active_scalar
        self.xy_pc = self.ax1.tricontourf(pts[:,0], pts[:,1], tri, val, **self.pc_props)

        # Draw Z-slice intersection in X-Z plot
        self.xz_ahw = self.ax2.axhline(self.zind, **self.clpropsw)
        self.xz_ahk = self.ax2.axhline(self.zind, **self.clpropsk)

        # Draw Z-slice intersection in Z-Y plot
        self.zy_avw = self.ax3.axvline(self.zind, **self.clpropsw)
        self.zy_avk = self.ax3.axvline(self.zind, **self.clpropsk)

    def update_xz(self):
        """Update plot for change in Y-index."""
        print("updating xz")
        # Clean up
        self._clear_elements(['xz_pc', 'zy_ahk', 'zy_ahw', 'xy_ahk', 'xy_ahw'])

        # Draw X-Z slice
        cntr = self.mesh.center
        cntr[1] = self.yind
        self.xz_slc = self.mesh.slice(normal="y", origin=cntr, generate_triangles=True)
        pts = self.xz_slc.points
        tri = self.xz_slc.faces.reshape((-1,4))[:, 1:]
        val = self.xz_slc.active_scalar
        self.xz_pc = self.ax2.tricontourf(pts[:,0], pts[:,2], tri, val, **self.pc_props)

        # Draw X-slice intersection in X-Y plot
        self.xy_ahw = self.ax1.axhline(self.yind, **self.clpropsw)
        self.xy_ahk = self.ax1.axhline(self.yind, **self.clpropsk)

        # Draw X-slice intersection in Z-Y plot
        self.zy_ahw = self.ax3.axhline(self.yind, **self.clpropsw)
        self.zy_ahk = self.ax3.axhline(self.yind, **self.clpropsk)

    def update_zy(self):
        """Update plot for change in X-index."""
        print("updating zy")
        # Clean up
        self._clear_elements(['zy_pc', 'xz_avw', 'xz_avk', 'xy_avw', 'xy_avk'])

        # Draw Z-Y slice
        cntr = self.mesh.center
        cntr[0] = self.xind
        self.zy_slc = self.mesh.slice(normal="x", origin=cntr, generate_triangles=True)
        pts = self.zy_slc.points
        tri = self.zy_slc.faces.reshape((-1,4))[:, 1:]
        val = self.zy_slc.active_scalar
        self.zy_pc = self.ax3.tricontourf(pts[:,2], pts[:,1], tri, val, **self.pc_props)

        # Draw Y-slice intersection in X-Y plot
        self.xy_avw = self.ax1.axvline(self.xind, **self.clpropsw)
        self.xy_avk = self.ax1.axvline(self.xind, **self.clpropsk)

        # Draw Y-slice intersection in X-Z plot
        self.xz_avw = self.ax2.axvline(self.xind, **self.clpropsw)
        self.xz_avk = self.ax2.axvline(self.xind, **self.clpropsk)

    def _clear_elements(self, names):
        """Remove elements from list <names> from plot if they exists."""
        for element in names:
            if hasattr(self, element):
                try:
                    getattr(self, element).remove()
                except:
                    pass


    @property
    def slices(self):
        return pv.MultiBlock([self.xy_slc, self.xz_slc, self.zy_slc])
