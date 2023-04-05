import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from uncertainties.unumpy import nominal_values as noms

def escatter(x, y, ax=None, **kwargs):
    """"
    Function for making a scatter plot with uncertainties.unumpy.uarray data.
    
    Parameters
    ----------
    x : uncerainty.unumpy.uarray or numpy.ndarray
    y : uncerainty.unumpy.uarray or numpy.ndarray
    ax : matplotlib.axes.Axes
    kwargs
        passed to plt.scatter
    
    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(1,1)
    if isinstance(x, pd.core.frame.DataFrame):
        x = x.values.flatten()
    if isinstance(y, pd.core.frame.DataFrame):
        y = y.values.flatten()
        
    if 'c' in kwargs:
        ma = ax.scatter(unp.nominal_values(x), unp.nominal_values(y), **kwargs)
    else:
        ma = None
        ax.scatter(unp.nominal_values(x), unp.nominal_values(y), **kwargs)
    
    ax.errorbar(unp.nominal_values(x), unp.nominal_values(y),
                xerr=unp.std_devs(x), 
                yerr=unp.std_devs(y),
                lw=0, elinewidth=1, color=(0,0,0,0.4), zorder=-1)
    
    return ax, ma

def pick_cax_position(x, y, f_width, f_height):
    x = noms(x)
    y = noms(y)
    xr = np.ptp(x)
    yr = np.ptp(y)
    
    xclear = xr * f_width
    yclear = yr * f_height
    
    # check number of points in each corner
    left = x < x.min() + xclear
    right = x > x.max() - xclear
    lower = y < y.min() + yclear
    upper = y > y.max() - yclear
    
    pos = np.array([
        np.sum(upper & right),  # 1
        np.sum(upper & left),  # 2
        np.sum(lower & left),  # 3
        np.sum(lower & right),  # 4
    ])
    
    # return choice in order of preference    
    preference = [4, 1, 3, 2]
    for p in preference:
        if pos[p-1] == 0:
            return p
    
    # otherwise just return the best one
    return np.argwhere(pos == min(pos)).flatten()[0] + 1
    
def unit_picker(a):
    base = np.log10(a)
    if base > -1:
        return 1, 'mol/mol'
    elif base > -4:
        return 1e3, 'mmol/mol'
    elif base > -7:
        return 1e6, 'µmol/mol'
    else:
        return 1e9, 'nmol/mol'
    
def match_lim(x, y, ax, pad=0.1, force_zero_min=False, line=True):
    x = noms(x)
    y = noms(y)
    
    if force_zero_min:
        lim = [np.max([0.0, np.min([x, y])]), np.max([x, y])]
    else:
        lim = [np.min([x, y]), np.max([x, y])]
    
    pad = np.array([-pad, pad]) * np.ptp(lim)
        
    lim += pad
        
    ax.set_xlim(lim)
    ax.set_ylim(lim)
        
    if line:
        ax.plot(lim, lim, ls='--', c='k', lw=1, alpha=0.4, zorder=-10)

expmap = {
    'Mg/Ca': ['Mg', 'Mg+Sr'],
    'Sr/Ca': ['Sr', 'Mg+Sr'],
    'B/C': ['B'],
    'B/Ca': ['B'],
    'Na/Ca': ['Mg', 'Sr', 'Mg+Sr', 'B']
    }

expcolors = {
    'Mg': 'C0',
    'Mg+Sr': 'C1',
    'Sr': 'C2',
    'B': 'C3',
}

varcmaps = {
    'Mg/Ca': plt.cm.Greens,
    'Sr/Ca': plt.cm.Blues,
}

ratio_to_partitioning = {
    'Mg/Ca': 'D_Mg',
    'Sr/Ca': 'D_Sr',
    'B/Ca': 'D_B',
    'B/C': 'D_B',
    'Na/Ca': 'D_Na',
}


def solution_vs_solid(dat, vars=['Mg/Ca', 'Sr/Ca', 'B/C', 'Na/Ca'], phase='overgrowth', xmode='solution_start', solid_mode='ratios', panel_size=3, axs=None):

    if axs is None:
        n = len(vars)
        fig, axs = plt.subplots(1, n, figsize=[panel_size * n, panel_size], constrained_layout=True)
    else:
        fig = axs[0].figure

    cind = (dat.metadata.Experiment.NA == 'Control').values.ravel()

    for var, ax in zip(vars, axs):
        
        if var == 'B/C':
            yvar = 'B/Ca'
        else:
            yvar = var
        
        if solid_mode != 'ratios':
            yvar = ratio_to_partitioning[yvar]
        
        # m, unit = unit_picker(np.quantile(noms(dat.loc[:, (phase, yvar)]), 0.5))
        # print(m, unit)
        m, unit = 1e3, 'mmol/mol'
        
        if var in ['Mg/Ca', 'Sr/Ca']:
            ind = np.zeros(dat.shape[0], dtype=bool)
            for exp in expmap[var]:
                ind = ind | (dat.metadata.Experiment.NA == exp).values.ravel()

            if var == 'Mg/Ca':
                cvar = 'Sr/Ca'
                c = noms(dat.loc[ind, (xmode, cvar)])
                vmax = np.nanmax(noms(dat.loc[:, (xmode, cvar)]))
            else:
                cvar = 'Mg/Ca'
                c = noms(dat.loc[ind, (xmode, cvar)])
                vmax = np.nanmax(noms(dat.loc[:, (xmode, cvar)]))
            
            cmap = varcmaps[var]
                    
            _, ma = escatter(
                x=dat.loc[ind, (xmode, var)].values[:,0], 
                y=dat.loc[ind, (phase, yvar)].values[:,0] * m,
                ax=ax, c=c, cmap=cmap, label=exp,
                vmin=0, vmax=vmax,
                edgecolor=(.4,.4,.4), lw=1)
            
            cpos = pick_cax_position(
                x=dat.loc[ind, (xmode, var)].values[:,0], 
                y=dat.loc[ind, (phase, yvar)].values[:,0] * m,
                f_width=0.3, f_height=0.2)
            cax = inset_axes(ax, width='30%', height='5%', loc=cpos)
            plt.colorbar(ma, cax=cax, label=f'Solution\n{cvar}', orientation='horizontal')
            if cpos in [3,4]:
                cax.xaxis.set_ticks_position('top')
                cax.xaxis.set_label_position('top')
            
        else:
            for exp in expmap[var]:
                ind = (dat.metadata.Experiment == exp).values            
                escatter(
                    x=dat.loc[ind, (xmode, var)].values[:,0], 
                    y=dat.loc[ind, (phase, yvar)].values[:,0] * m,
                    ax=ax, c=expcolors[exp], label=exp,
                    edgecolor=(.4,.4,.4), lw=1)
                
        escatter(
            x=dat.loc[cind, (xmode, var)].values[:,0], 
            y=dat.loc[cind, (phase, yvar)].values[:,0] * m,
            ax=ax, c='k', marker='s', label='Control',
            edgecolor=(.4,.4,.4), lw=1)
        
        ax.set_xlabel(f'{var} solution (mol/mol)')
        ax.set_ylabel(f'{yvar} {phase} ({unit})')

    axs[-1].legend(title='Experiment')

    return fig, axs

def solid_vs_solid(dat, vars=['Mg/Ca', 'Sr/Ca', 'B/C', 'Na/Ca'], xphase='overgrowth', yphase='overgrowth', cmode='solution_start', solid_mode='ratios', panel_size=3, match_axes=True, axs=None):

    if axs is None:
        n = len(vars)
        fig, axs = plt.subplots(1, n, figsize=[panel_size * n, panel_size], constrained_layout=True)
    else:
        fig = axs[0].figure

    cind = (dat.metadata.Experiment.NA == 'Control').values.ravel()

    for var, ax in zip(vars, axs):
        
        if solid_mode != 'ratios':
            var = ratio_to_partitioning[var]
        
        xm, xunit = unit_picker(np.nanmean(noms(dat.loc[:, (xphase, var)])))
        m, unit = unit_picker(np.nanmean(noms(dat.loc[:, (yphase, var)])))
        
        if xm < m:
            m = xm
            unit = xunit
        
        if var in ['Mg/Ca', 'Sr/Ca']:
            ind = np.zeros(dat.shape[0], dtype=bool)
            for exp in expmap[var]:
                ind = ind | (dat.metadata.Experiment.NA == exp).values.ravel()

            if var == 'Mg/Ca':
                cvar = 'Sr/Ca'
                c = noms(dat.loc[ind, (cmode, cvar)])
                vmax = np.nanmax(noms(dat.loc[:, (cmode, cvar)]))
            else:
                cvar = 'Mg/Ca'
                c = noms(dat.loc[ind, (cmode, cvar)])
                vmax = np.nanmax(noms(dat.loc[:, (cmode, cvar)]))
            
            cmap = varcmaps[var]
                    
            _, ma = escatter(
                x=dat.loc[ind, (xphase, var)].values[:,0] * m, 
                y=dat.loc[ind, (yphase, var)].values[:,0] * m,
                ax=ax, c=c, cmap=cmap, label=exp,
                vmin=0, vmax=vmax,
                edgecolor=(.4,.4,.4), lw=1)
            
            cpos = pick_cax_position(
                x=dat.loc[ind, (xphase, var)].values[:,0] * m, 
                y=dat.loc[ind, (yphase, var)].values[:,0] * m,
                f_width=0.3, f_height=0.2)
            cax = inset_axes(ax, width='30%', height='5%', loc=cpos)
            plt.colorbar(ma, cax=cax, label=f'Solution\n{cvar}', orientation='horizontal')
            if cpos in [3,4]:
                cax.xaxis.set_ticks_position('top')
                cax.xaxis.set_label_position('top')
            
            if match_axes:
                match_lim(x=dat.loc[ind, (xphase, var)].values[:,0] * m, 
                          y=dat.loc[ind, (yphase, var)].values[:,0] * m,
                          ax=ax)
            
        else:
            for exp in expmap[var]:
                ind = (dat.metadata.Experiment == exp).values            
                escatter(
                    x=dat.loc[ind, (xphase, var)].values[:,0] * m, 
                    y=dat.loc[ind, (yphase, var)].values[:,0] * m,
                    ax=ax, c=expcolors[exp], label=exp,
                    edgecolor=(.4,.4,.4), lw=1)
                
            if match_axes:
                match_lim(x=dat.loc[ind, (xphase, var)].values[:,0] * m, 
                          y=dat.loc[ind, (yphase, var)].values[:,0] * m,
                          ax=ax)
        
        escatter(
            x=dat.loc[cind, (xphase, var)].values[:,0] * m, 
            y=dat.loc[cind, (yphase, var)].values[:,0] * m,
            ax=ax, c='k', marker='s', label='Control',
            edgecolor=(.4,.4,.4), lw=1)
        
        ax.set_xlabel(f'{var} {xphase} ({unit})')
        ax.set_ylabel(f'{var} {yphase} ({unit})')

    axs[-1].legend(title='Experiment')
    
    return fig, axs

def solution_vs_yvar(dat, vars=['Mg/Ca', 'Sr/Ca', 'B/C', 'Na/Ca'], yvar=('overgrowth', 'F_V'), xmode='solution_start', panel_size=3, axs=None):

    if axs is None:
        n = len(vars)
        fig, axs = plt.subplots(1, n, figsize=[panel_size * n, panel_size], constrained_layout=True)
    else:
        fig = axs[0].figure
    
    cind = (dat.metadata.Experiment.NA == 'Control').values.ravel()
    
    for var, ax in zip(vars, axs):
                    
        if var in ['Mg/Ca', 'Sr/Ca']:
            ind = np.zeros(dat.shape[0], dtype=bool)
            
            for exp in expmap[var]:
                ind = ind | (dat.metadata.Experiment.NA == exp).values.ravel()

            if var == 'Mg/Ca':
                cvar = 'Sr/Ca'
                c = noms(dat.loc[ind, (xmode, cvar)])
                vmax = np.nanmax(noms(dat.loc[:, (xmode, cvar)]))
            else:
                cvar = 'Mg/Ca'
                c = noms(dat.loc[ind, (xmode, cvar)])
                vmax = np.nanmax(noms(dat.loc[:, (xmode, cvar)]))
            
            cmap = varcmaps[var]
                    
            _, ma = escatter(
                x=dat.loc[ind, (xmode, var)].values[:,0], 
                y=dat.loc[ind, yvar].values[:,0],
                ax=ax, c=c, cmap=cmap, label=exp,
                vmin=0, vmax=vmax,
                edgecolor=(.4,.4,.4), lw=1)
            
            cpos = pick_cax_position(
                x=dat.loc[ind, (xmode, var)].values[:,0], 
                y=dat.loc[ind, yvar].values[:,0],
                f_width=0.3, f_height=0.25)
            cax = inset_axes(ax, width='30%', height='5%', loc=cpos)
            plt.colorbar(ma, cax=cax, label=f'Solution\n{cvar}', orientation='horizontal')
            if cpos in [3,4]:
                cax.xaxis.set_ticks_position('top')
                cax.xaxis.set_label_position('top')
            
        else:
            for exp in expmap[var]:
                ind = (dat.metadata.Experiment == exp).values            
                escatter(
                    x=dat.loc[ind, (xmode, var)].values[:,0], 
                    y=dat.loc[ind, yvar].values[:,0],
                    ax=ax, c=expcolors[exp], label=exp,
                    edgecolor=(.4,.4,.4), lw=1)
        
        escatter(
            x=dat.loc[cind, (xmode, var)].values[:,0], 
            y=dat.loc[cind, yvar].values[:,0],
            ax=ax, c='k', marker='s',
            edgecolor=(.4,.4,.4), lw=1)
        
        ax.set_xlabel(f'{var} solution (mol/mol)')
        ax.set_ylabel(f'{yvar[0]} {yvar[1]}')

    axs[-1].legend()

    return fig, axs

def label_axes(axs, labels=None, x=0.02, y=0.98, va='top', ha='left', fontsize=12, weight='bold', color='grey', zorder=999, **kwargs):
    if labels is None:
        labels = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for ax, lab in zip(axs.flat, labels):
        ax.text(x, y, lab, transform=ax.transAxes, va=va, ha=ha, fontsize=fontsize, weight=weight, color=color, zorder=zorder, **kwargs)
        