#!/usr/bin/python
# _*_ coding: latin-1 -*-

import numpy as np
from numpy import ma
import mpl_toolkits.axisartist.grid_finder as GF
import mpl_toolkits.axisartist.floating_axes as FA
import matplotlib.pyplot as plt
import netCDF4

from matplotlib.projections.polar import PolarAxes

def make_Taylor_axis(smax=1.5, figsize=None, fig=None):
    # radial grid
    # rlocs = np.concatenate((np.arange(0,1,0.2),[0.95,0.99]))  # correlation coefficient tick marks
    rlocs = [0, .3, .6, .8, .95, .99]
    str_rlocs = np.concatenate((np.arange(0,1,0.2), [0.95,0.99], np.arange(0,10,0.2),[0.95,0.99]))
    tlocs = np.arccos(rlocs)        # Conversion to polar angles
    gl1 = GF.FixedLocator(tlocs)    # Positions
    tf1 = GF.DictFormatter(dict(zip(tlocs, map("{:.2f}".format, rlocs))))

    # circumferencial grid
    str_locs2 = tlocs2 =  np.arange(0, 2, 0.2)
    g22 = GF.FixedLocator(tlocs2)  
    tf2 = GF.DictFormatter(dict(zip(tlocs2, map("{:.1f}".format, str_locs2))))

    tr = PolarAxes.PolarTransform()
    
    smin = 0

    ghelper = FA.GridHelperCurveLinear(tr,
                                       extremes=(0, np.pi/2, # 1st quadrant
                                                 smin,smax),
                                       grid_locator1=gl1,
                                       #grid_locator2=g11,
                                       tick_formatter1=tf1,
                                       tick_formatter2=tf2,
                                       )
    if fig is None:
        fig = plt.figure(figsize=figsize, dpi=100)
    
    ax = FA.FloatingSubplot(fig, 121, grid_helper=ghelper)

    fig.add_subplot(ax)
    ax.axis["top"].set_axis_direction("bottom") 
    ax.axis["top"].toggle(ticklabels=True, label=True)
    ax.axis["top"].major_ticklabels.set_axis_direction("top")
    ax.axis["top"].label.set_axis_direction("top")
    ax.axis["top"].label.set_text("Correlation Coefficient")

    ax.axis["left"].set_axis_direction("bottom") 
    ax.axis["left"].label.set_text("Standard Deviation")

    ax.axis["right"].set_axis_direction("top") 
    ax.axis["right"].toggle(ticklabels=True, label=True)
    ax.axis["right"].set_visible(True)
    ax.axis["right"].major_ticklabels.set_axis_direction("bottom")
    #ax.axis["right"].label.set_text("Standard Deviation")

    ax.axis["bottom"].set_visible(False) 

    ax.grid(True)

    ax = ax.get_aux_axes(tr)

    # radial RMS error lines from observed.
    ref = 1
    t = np.linspace(0, np.pi)
    r = np.zeros_like(t) + ref
    ax.plot(t, r, 'k--', label='_')

    rs,ts = np.meshgrid(np.linspace(smin, smax),
                        np.linspace(0,np.pi))
    
    rms = np.sqrt(ref**2 + rs**2 - 2*ref*rs*np.cos(ts))
    CS = ax.contour(ts, rs, rms, colors='k', alpha=0.6, linewidths=1)
    plt.clabel(CS, inline=1, fontsize=10)
    
    # if figsize is None:
    #     fig.set_size_inches(5, 5)
    # else:
    #     fig.set_size_inches(*figsize)
    
    return fig, ax

def Taylor_diag(series, ax=None, st_dev_max=None,  **kwargs):
    """ Taylor Diagram : obs is reference data sample
        in a full diagram (0 --> npi)
        --------------------------------------------------------------------------
        Input: series     - list with all time series (arrays) to analyze  
               series[0]  - is the observation, the reference by default.
    """
    corr = np.zeros(len(series))
    std = np.zeros(len(series))

    for i in range(len(series)):
        corr[i] = ma.corrcoef(series[0],series[i])[1,0]
        std[i] = ma.std(series[i])/ma.std(series[0])

    ref = 1# ma.std(series[0])
    
    if ax is None:
        _, ax = make_Taylor_axis()
    
    # plot observed point
    ax.plot(np.arccos(0.9999),ref,'k',marker='*',ls='', ms=10)

    # plot modelled points
    if st_dev_max is not None:
        ind = std < st_dev_max
        ind[0] = True
        corr = corr[ind]
        std = std[ind]
        if "c" in kwargs:
            kwargs["c"] = kwargs["c"][ind[1:]]

    point = ax.scatter(np.arccos(corr[1:]), std[1:], **kwargs)
        # ax.text(np.arccos(corr[i]), std[i], i)
    
    # plt.legend(bbox_to_anchor=(1.5, 1),prop=dict(size='large'),loc='best', fontsize=9)
    
    return point

