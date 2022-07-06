
import numpy as np
import pandas as pd
from scipy.special import voigt_profile
from scipy.optimize import curve_fit
import uncertainties as un

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

def fit_peak(dat, location, window=0.5):

    sub = dat.loc[(dat.Angle >= location - window) & (dat.Angle <= location + window)]
    p0 = (
        sub.loc[sub.PSD == max(sub.PSD), 'Angle'].item(), 
        sub.loc[abs(sub.Angle - location) == min(abs(sub.Angle - location)), 'PSD'].item() - sub['PSD'].iloc[0],
        .1,
        .1,
        sub['PSD'].iloc[0],
        0
    )
    
    return curve_fit(pkfn, sub.Angle, sub.PSD, p0=p0, maxfev=2000)

def calc_vat_frac(file, vat_location=12.42, cal_location=13.47, window=0.5):
    
    dat = read_xrd(file)
        
    vat_peak = fit_peak(dat, location=vat_location, window=window)
    cal_peak = fit_peak(dat, location=cal_location, window=window)
    
    areas = []
    for p, cov in [vat_peak, cal_peak]:
        xn = np.linspace(p[0] - window, p[0] + window)
        yn = pkfn(xn, *p)
        yn -= p[-2]
        yn -= p[-1] * (xn - p[0])
        areas.append(np.trapz(yn, xn))
    
    R = areas[0] / areas[1] / OFFSET
    
    return R / (1 + R)