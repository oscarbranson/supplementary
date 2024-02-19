import os
import numpy as np
import uncertainties as un
from uncertainties.unumpy import nominal_values, std_devs
from scipy.optimize import curve_fit

rdict = {
    'This Study': 'This Study',
    'Gabitov2019': 'Gabitov et al. (2019)',
    'Langer2015': 'Langer et al. (2015)', 
    'Marriott2004a': 'Marriott et al. (2004a)', 
    'Marriott2004b': 'Marriott et al. (2004b)', 
    'Roberts2018': 'Roberts et al. (2018)', 
    'Vigier2015': 'Vigier et al. (2015)',
    'Fuger2019': 'Fuger et al. (2019, 2022)',
    'Day2021': 'Day et al. (2021)',
}

mdict = {
    'This Study': ('C0', 'o'),
    'Gabitov2019': ('C2', '^'),
    'Langer2015': ('C5', 'd'),
    'Marriott2004a': ('C4', '<'),
    'Marriott2004b': ('C4', '>'),
    'Roberts2018': ('C6', '+'),
    'Vigier2015': ('C7', 'x'),
    'Fuger2019': ('C3', 's'),
    'Day2021': ('C2', 'v'),
}

def savefig(fig, name, figdir='figs'):
    fig.savefig(os.path.join(figdir, name + '.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(figdir, name + '.pdf'))
    
def line(x, m, c):
    return x * m + c

def fit_line(x, y):
    nona = ~np.isnan(x) & ~np.isnan(y)
    
    return un.correlated_values(*curve_fit(line, x[nona], y[nona]))

def plot_line(ax, coef, xn=None, **kwargs):
    xlim = ax.get_xlim()
    ax.set_xlim(xlim)
    
    if xn is None:
        xn = np.linspace(*xlim, 50)
    yn = line(xn, *coef)
    
    ax.plot(xn, nominal_values(yn), **kwargs)
    if 'ls' in kwargs:
        kwargs.pop('ls')
    ax.fill_between(xn, nominal_values(yn) - 1.96 * std_devs(yn), nominal_values(yn) + 1.96 * std_devs(yn), alpha=0.2, **kwargs, lw=0)