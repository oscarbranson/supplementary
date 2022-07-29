
from ctypes.wintypes import VARIANT_BOOL
import numpy as np
import pandas as pd
from pyproj import transform
from scipy.special import voigt_profile
from scipy.optimize import curve_fit
import uncertainties as un
import matplotlib.pyplot as plt

# parameter calculated in XRD_Calibration
OFFSET = un.ufloat(0.4654043388930521, 0.007372891362069755)


def read_xrd(file):
    with open(file, 'r') as f:
        i = 0
        foundstart = False
        while not foundstart:
            if '[Data]' in f.readline():  # [Data] marks the start of the data
                foundstart = True
            i += 1
    
    dat = pd.read_csv(file, skiprows=i).dropna(axis=1)
    dat.columns = dat.columns.str.strip()
    return dat

def pkfn(angle, position, height, sigma, gamma, background, background_slope):
    shape = voigt_profile(angle - position, sigma, gamma)
    shape /= shape.max()
    return height * shape + background + (angle - position) * background_slope

def pkfn_nobkg(angle, position, height, sigma, gamma, background, background_slope):
    shape = voigt_profile(angle - position, sigma, gamma)
    shape /= shape.max()
    return height * shape

def fit_peak(dat, location, win_left=0.5, win_right=0.5, update_loc=True):
    
    sub = dat.loc[(dat.Angle >= location - win_left) & (dat.Angle <= location + win_right)]
    
    # find max
    if update_loc:
        location = (location + sub.loc[sub.PSD == max(sub.PSD), 'Angle'].item()) / 2
        sub = dat.loc[(dat.Angle >= location - win_left) & (dat.Angle <= location + win_right)]
    
    p0 = (
        sub.loc[sub.PSD == max(sub.PSD), 'Angle'].item(), 
        sub.loc[abs(sub.Angle - location) == min(abs(sub.Angle - location)), 'PSD'].item() - sub['PSD'].iloc[0],
        .1,
        .1,
        sub['PSD'].iloc[0],
        0
    )
    
    return curve_fit(pkfn, sub.Angle, sub.PSD, p0=p0, maxfev=2000)

def calc_vat_frac(file, vat_location=12.42, cal_location=13.47, win_left=0.5, win_right=None, a=OFFSET):
    
    if win_right is None:
        win_right = win_left 
        
    dat = read_xrd(file)
        
    vat_peak = fit_peak(dat, location=vat_location, win_left=win_left, win_right=win_right)
    cal_peak = fit_peak(dat, location=cal_location, win_left=win_left, win_right=win_right)
    
    areas = []
    for p, cov in [vat_peak, cal_peak]:
        xn = np.linspace(p[0] - win_left, p[0] + win_right)
        yn = pkfn(xn, *p)
        yn -= p[-2]
        yn -= p[-1] * (xn - p[0])
        areas.append(np.trapz(yn, xn))
    
    R = areas[0] / areas[1] / a
    
    return R / (1 + R)

def plot_fit(file, vat_location=12.42, cal_location=13.47, win_left=0.5, win_right=None, update_loc=True):
    dat = read_xrd(file)
    
    if win_right is None:
        win_right = win_left 
        
    vat_peak = fit_peak(dat, location=vat_location, win_left=win_left, win_right=win_right, update_loc=update_loc)
    cal_peak = fit_peak(dat, location=cal_location, win_left=win_left, win_right=win_right, update_loc=update_loc)
    
    if update_loc:
        vat_location = vat_peak[0][0]
        cal_location = cal_peak[0][0]
    
    fig, ax = plt.subplots(figsize=(6, 2.4), constrained_layout=True)
    ax.scatter(dat.Angle, dat.PSD, s=5, color=(0,0,0,0.5))
    ax.set_xlim(11.5, 14.5)
    
    for location, peak, c in zip([vat_location, cal_location], [vat_peak, cal_peak], ['C1', 'C2']):
    
        ax.axvline(location, color=c)
        ax.axvspan(location - win_left, location + win_right, color=c, alpha=0.2)
        
        xn = np.linspace(location - win_left, location + win_right)
        ax.plot(xn, pkfn(xn, *peak[0]), color=c)
    
    peak = vat_peak[0]
    ax.text(0.02, 0.98, f'loc: {peak[0]:.2f}\nheight: {peak[1]:.2f}\nsigma: {peak[2]:.2f}\ngamma: {peak[3]:.2f}\nbkg_c: {peak[4]:.2f}\nbkg_m: {peak[5]:.2f}', va='top', ha='left', transform=ax.transAxes, color='C1', fontsize=8)

    peak = cal_peak[0]
    ax.text(0.98, 0.98, f'loc: {peak[0]:.2f}\nheight: {peak[1]:.2f}\nsigma: {peak[2]:.2f}\ngamma: {peak[3]:.2f}\nbkg_c: {peak[4]:.2f}\nbkg_m: {peak[5]:.2f}', va='top', ha='right', transform=ax.transAxes, color='C2', fontsize=8)
    
    return fig, ax