import numpy as np
from scipy.optimize import curve_fit
from uncertainties import unumpy as unp

from .helpers import read_data, isolate_constant_conditions

nom = unp.nominal_values
err = unp.std_devs

# previously proposed relationship functions
### Allen et al 2011 (Orbulina)
def bca_Allen2011(BO4, m=1.4, c=39.5):
    """all units in umol"""
    return m * BO4 + c

### Foster et al 2008 (ruber)
def kb_Foster2008_CO3(CO3, m=-0.00842, c=3.9708):
    """input in umol"""
    return (m * CO3 + c) / 1000

### Henehan 2015
def bca_Henehan2015_pH(pH, m=104.64, c=667.56):
    return pH * m + c

def bca_Henehan2015_BO4_HCO3(BO4_HCO3, m=1104.34, c=102.76):
    return BO4_HCO3 * m + c

def bca_Henehan2015_BO4_DIC(BO4_DIC, m=1378.2, c=97.51):
    return BO4_DIC * m + c

# construct data for fitting
def get_training_subset(dat):
    return isolate_constant_conditions(dat, 
                                  {
                                      ('csys_mid', 'DIC'): (2050, 100),
                                      ('Measured', '[Mg]sw'): (50, 5),
                                      ('Measured', '[Ca]sw'): (10, 1),
                                  })

def make_fitting_data(dat):
    functions = {}
    xdata_all = {}
    xdata_training_conditions = {}
    ydata_all = {}
    ydata_training_conditions = {}

    sub = get_training_subset(dat)

    functions['Allen_2011_bca'] = bca_Allen2011
    xdata_all['Allen_2011_bca'] = dat.csys_mid.BO4
    xdata_training_conditions['Allen_2011_bca'] = sub.csys_mid.BO4
    ydata_all['Allen_2011_bca'] = nom(dat.loc[:, ('Measured', 'B/Caf')])
    ydata_training_conditions['Allen_2011_bca'] = nom(sub.loc[:, ('Measured', 'B/Caf')])

    functions['Foster2008_kb'] = kb_Foster2008_CO3
    xdata_all['Foster2008_kb'] = dat.csys_mid.CO3
    xdata_training_conditions['Foster2008_kb'] = sub.csys_mid.CO3
    ydata_all['Foster2008_kb'] = nom(dat.loc[:, ('Measured', 'KB')])
    ydata_training_conditions['Foster2008_kb'] = nom(sub.loc[:, ('Measured', 'KB')])

    functions['Henehan2015_pH_bca'] = bca_Henehan2015_pH
    xdata_all['Henehan2015_pH_bca'] = dat.csys_mid.pHtot
    xdata_training_conditions['Henehan2015_pH_bca'] = sub.csys_mid.pHtot
    ydata_all['Henehan2015_pH_bca'] = nom(dat.loc[:, ('Measured', 'B/Caf')])
    ydata_training_conditions['Henehan2015_pH_bca'] = nom(sub.loc[:, ('Measured', 'B/Caf')])

    functions['Henehan2015_BO4_HCO3_bca'] = bca_Henehan2015_BO4_HCO3
    xdata_all['Henehan2015_BO4_HCO3_bca'] = dat.csys_mid.BO4 / dat.csys_mid.CO3
    xdata_training_conditions['Henehan2015_BO4_HCO3_bca'] = sub.csys_mid.BO4 / sub.csys_mid.CO3
    ydata_all['Henehan2015_BO4_HCO3_bca'] = nom(dat.loc[:, ('Measured', 'B/Caf')])
    ydata_training_conditions['Henehan2015_BO4_HCO3_bca'] = nom(sub.loc[:, ('Measured', 'B/Caf')])

    functions['Henehan2015_BO4_DIC_bca'] = bca_Henehan2015_BO4_DIC
    xdata_all['Henehan2015_BO4_DIC_bca'] = dat.csys_mid.BO4 / dat.csys_mid.DIC
    xdata_training_conditions['Henehan2015_BO4_DIC_bca'] = sub.csys_mid.BO4 / sub.csys_mid.DIC
    ydata_all['Henehan2015_BO4_DIC_bca'] = nom(dat.loc[:, ('Measured', 'B/Caf')])
    ydata_training_conditions['Henehan2015_BO4_DIC_bca'] = nom(sub.loc[:, ('Measured', 'B/Caf')])

    return functions, xdata_all, ydata_all, xdata_training_conditions, ydata_training_conditions

def fit_data_with_previous_models(dat):
    fits_all_data = {}
    fits_training_conditions = {}

    preds_all_data = {}
    preds_training_conditions = {}

    functions, xdata_all, ydata_all, xdata_training_conditions, ydata_training_conditions = make_fitting_data(dat)

    for k, f in functions.items():
        p, _ = curve_fit(f, xdata_all[k], ydata_all[k])
        fits_all_data[k] = p
        preds_all_data[k] = f(xdata_all[k], *p)
        
        p_train, _ = curve_fit(f, xdata_training_conditions[k], ydata_training_conditions[k])
        fits_training_conditions[k] = p_train
        preds_training_conditions[k] = f(xdata_training_conditions[k], *p_train)
    
    return {
        "fits_all_data": fits_all_data, 
        "preds_all_data": preds_all_data, 
        "ydata_all": ydata_all,
        "fits_training_conditions": fits_training_conditions, 
        "preds_training_conditions": preds_training_conditions,
        "ydata_training_conditions": ydata_training_conditions
        }