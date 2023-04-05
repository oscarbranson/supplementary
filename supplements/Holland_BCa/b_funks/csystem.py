"""
Functions for calculating the C system using pitzer and MyAMI_V1
"""
import pandas as pd
import cbsyst as cb  # must be co2sys_compare branch
from tqdm import tqdm
from .phreeq import calc_pitzer_fKs

def calc_pitzer_Csys(d, 
                     database='./b_funks/db_phreeqc/pitzer.dat', 
                     phreeq_path='/usr/local/lib/libiphreeqc.so'):
    
    # define standard seawater composition (mol/kg SW)
    # ================================================
    sw = {'pH': 8.2,
          'Temp': 25,
          'units': 'mol/kgw',
          'density': 1.026,
          'Ca': 10.28e-3,
          'Cl': 0.54586,
          'K': 10.21e-3,
          'Mg': 52.82e-3,
          'Na': 0.46906,
          'S(6)': 0.02824,
          'B': 416e-6,
          'C': 2.3e-3,
          'Br': 0.84e-3,
          'F': 0.7e-4}
    
    # calculate ambient Ks at 25C (measurement temperature)
    # =====================================================
    # ks = {}

    TempC = 25
    TempK = TempC + 273.15
    Sal = d.loc[:, ('Measured', 'Salinity')]

    ks = cb.calc_Ks(temp_c=TempC, sal=Sal)
    # par = cb.MyAMI_V2.start_params
    # fns = cb.MyAMI_V2.fn_dict
    
    # for k in ['K0', 'K1', 'K2', 'KB', 'KW', 'KspC', 'KspA', 'KSO4']:
    #     ks[k] = fns[k]((TempK, Sal), *par[k])

    # calculate other Ks
    # ks.update(cb.non_MyAMI_constants.calc_KF(TempC, Sal))
    # ks.update(cb.non_MyAMI_constants.calc_KPs(TempC, Sal))
    # ks.update(cb.non_MyAMI_constants.calc_KSi(TempC, Sal))
    
    # calculate correction factors at 25C
    # ===================================
    fKs = {}
    isw = sw.copy()

    salsc = ['Cl', 'K', 'Na', 'S(6)', 'Br', 'F']

    for i, v in tqdm(d.iterrows(), total=d.shape[0],
                     desc='Calculating 25C Corrections'):
        v = v.Measured
        isw.update({'Mg': v.loc['[Mg]sw'] * 1e-3,
                    'Ca': v.loc['[Ca]sw'] * 1e-3,
                    'Temp': 25,
                    'B': v.loc['[B]sw'] * 1e-6})

        sal_f = v.Salinity / 35.
        isw.update({k: sw[k] * sal_f for k in salsc})

        fks = calc_pitzer_fKs(sw, isw, database=database, phreeq_path=phreeq_path)

        for k, f in fks.items():
            if k not in fKs:
                fKs[k] = []
            fKs[k].append(f)
    
    # apply correction factors to 25C Ks
    # ==================================
    cKs = ks.copy()
    for k, f in fKs.items():
        cKs[k] = ks[k] * f
    
    # convert all pH scales to Total
    # ==============================
    pHs = cb.calc_pH_scales(pHtot=None, pHfree=None, pHsws=None,
                            pHNBS=d.loc[:, ('Measured', 'pHNBS')],
                            ST=cb.calc_ST(Sal), FT=cb.calc_FT(Sal),
                            TempK=TempK, Sal=Sal, Ks=cb.Bunch(cKs))

    d.loc[:, ('pitzer_25C', 'pHTOTAL')] = d.loc[:, ('Measured', 'pHTOTAL')]

    ind = d.loc[:, ('Measured', 'pHTOTAL')].isnull()
    d.loc[ind, ('pitzer_25C', 'pHTOTAL')] = pHs['pHtot'][ind]
    
    # Calculate DIC/Alk at MEASURED conditions (25C)
    # ==============================================
    # calculate carbon system using 25C Ks, pH and either DIC or Alk
    carb = cb.Csys(d.loc[:, ('pitzer_25C', 'pHTOTAL')], TA=d.loc[:, ('Measured', 'Alk')], 
                   BT=d.loc[:, ('Measured', '[B]sw')], Ks=cKs)
    d.loc[:, ('pitzer_25C', 'DIC')] = carb.DIC

    carb = cb.Csys(d.loc[:, ('pitzer_25C', 'pHTOTAL')], DIC=d.loc[:, ('Measured', 'DIC')], 
                   BT=d.loc[:, ('Measured', '[B]sw')], Ks=cKs)
    d.loc[:, ('pitzer_25C', 'Alk')] = carb.TA

    # whichever of DIC or Alk wasn't calculated, transfer the
    # measured value to the pitzer_25C column for in-situ condition calculation
    dicnull = d.loc[:, ('pitzer_25C', 'DIC')].isnull()
    d.loc[dicnull, ('pitzer_25C', 'DIC')] = d.loc[dicnull, ('Measured', 'DIC')]

    alknull = d.loc[:, ('pitzer_25C', 'Alk')].isnull()
    d.loc[alknull, ('pitzer_25C', 'Alk')] = d.loc[alknull, ('Measured', 'Alk')]
    
    
    # calculate empirical Ks for ambient Mg and Ca, and experiment Sal and Temp
    # =========================================================================
    eKs = cb.calc_Ks(d.loc[:, ('Measured', 'Temp')], d.loc[:, ('Measured', 'Salinity')], 0, sw['Mg'], sw['Ca'],
                     cb.calc_ST(d.loc[:, ('Measured', 'Salinity')]), cb.calc_FT(d.loc[:, ('Measured', 'Salinity')]))
    
    # calculate correction factors at experimental temperature
    # ========================================================
    fKs = {}
    isw = sw.copy()

    salsc = ['Cl', 'K', 'Na', 'S(6)', 'Br', 'F']

    for i, v in tqdm(d.iterrows(), total=d.shape[0],
                     desc='Calculating Experimental T Corrections'):
        v = v.Measured
        isw.update({'Mg': v.loc['[Mg]sw'] * 1e-3,
                    'Ca': v.loc['[Ca]sw'] * 1e-3,
                    'Temp': v.Temp,
                    'B': v.loc['[B]sw'] * 1e-6})

        sal_f = v.Salinity / 35.
        isw.update({k: sw[k] * sal_f for k in salsc})

        fks = calc_pitzer_fKs(sw, isw, database=database, phreeq_path=phreeq_path)

        for k, f in fks.items():
            if k not in fKs:
                fKs[k] = []
            fKs[k].append(f)
    
    
    # apply correction factors
    for k, f in fKs.items():
        eKs[k] = ks[k] * f

    # calculate carbon system at experimental conditions using pitzer-calculated DIC/Alk pairs
    carb = cb.CBsys(DIC=d.loc[:, ('pitzer_25C', 'DIC')], TA=d.loc[:, ('pitzer_25C', 'Alk')],
                    BT=d.loc[:, ('Measured', '[B]sw')], Ks=eKs)

    # store outputs!
    for c in ['DIC', 'TA', 'CO3', 'HCO3', 'CO2', 'pHtot', 'BO4', 'BO3', 'BT']:
        d.loc[:, ('pitzer', c)] = carb[c]
    
    # return carb
    d.loc[:, ('pitzer', 'KspC')] = carb['Ks']['KspC']
    
    return d

def calc_MyAMI_Csys(d):
    out_cols =  ['DIC', 'TA', 'CO3', 'HCO3', 'CO2', 'pHtot', 'BO4', 'BO3', 'BT']
    
    # for pHNBS and Alk
    ind = ~d.loc[:, ('Measured', 'pHNBS')].isnull() & ~d.loc[:, ('Measured', 'Alk')].isnull()
    if any(ind):
        cs = cb.CBsys(pHNBS=d.loc[ind, ('Measured', 'pHNBS')], TA=d.loc[ind, ('Measured', 'Alk')],
                     T_in=25, T_out=d.loc[ind, ('Measured', 'Temp')], S_in=d.loc[ind, ('Measured', 'Salinity')],
                     Mg=d.loc[ind, ('Measured', '[Mg]sw')] / 1e3, Ca=d.loc[ind, ('Measured', '[Ca]sw')] / 1e3,
                     BT=d.loc[ind, ('Measured', '[B]sw')])
        for c in out_cols:
            d.loc[ind, ('MyAMI', c)] = cs[c]
        d.loc[ind, ('MyAMI', 'KspC')] = cs['Ks']['KspC']
    
    # for pHNBS and DIC
    ind = ~d.loc[:, ('Measured', 'pHNBS')].isnull() & ~d.loc[:, ('Measured', 'DIC')].isnull()
    if any(ind):
        cs = cb.CBsys(pHNBS=d.loc[ind, ('Measured', 'pHNBS')], DIC=d.loc[ind, ('Measured', 'DIC')],
                    T_in=25, T_out=d.loc[ind, ('Measured', 'Temp')], S_in=d.loc[ind, ('Measured', 'Salinity')],
                    Mg=d.loc[ind, ('Measured', '[Mg]sw')] / 1e3, Ca=d.loc[ind, ('Measured', '[Ca]sw')] / 1e3,
                    BT=d.loc[ind, ('Measured', '[B]sw')])
        for c in out_cols:
            d.loc[ind, ('MyAMI', c)] = cs[c]
        d.loc[ind, ('MyAMI', 'KspC')] = cs['Ks']['KspC']
    
    # for pHTOTAL and Alk
    ind = ~d.loc[:, ('Measured', 'pHTOTAL')].isnull() & ~d.loc[:, ('Measured', 'Alk')].isnull()
    if any(ind):
        cs = cb.CBsys(pHtot=d.loc[ind, ('Measured', 'pHTOTAL')], TA=d.loc[ind, ('Measured', 'Alk')],
                    T_in=25, T_out=d.loc[ind, ('Measured', 'Temp')], S_in=d.loc[ind, ('Measured', 'Salinity')],
                    Mg=d.loc[ind, ('Measured', '[Mg]sw')] / 1e3, Ca=d.loc[ind, ('Measured', '[Ca]sw')] / 1e3,
                    BT=d.loc[ind, ('Measured', '[B]sw')])
        for c in out_cols:
            d.loc[ind, ('MyAMI', c)] = cs[c]
        d.loc[ind, ('MyAMI', 'KspC')] = cs['Ks']['KspC']
    
    # for pHTOTAL and DIC
    ind = ~d.loc[:, ('Measured', 'pHTOTAL')].isnull() & ~d.loc[:, ('Measured', 'DIC')].isnull()
    if any(ind):
        cs = cb.CBsys(pHtot=d.loc[ind, ('Measured', 'pHTOTAL')], DIC=d.loc[ind, ('Measured', 'DIC')],
                    T_in=25, T_out=d.loc[ind, ('Measured', 'Temp')], S_in=d.loc[ind, ('Measured', 'Salinity')],
                    Mg=d.loc[ind, ('Measured', '[Mg]sw')] / 1e3, Ca=d.loc[ind, ('Measured', '[Ca]sw')] / 1e3,
                    BT=d.loc[ind, ('Measured', '[B]sw')])
        for c in out_cols:
            d.loc[ind, ('MyAMI', c)] = cs[c]
        d.loc[ind, ('MyAMI', 'KspC')] = cs['Ks']['KspC']
    
    return d

def mean_Csys(d):
    out_cols =  ['DIC', 'TA', 'CO3', 'HCO3', 'CO2', 'pHtot', 'BO4', 'BO3', 'BT', 'KspC']

    for c in out_cols:
        d.loc[:, ('csys_mid', c)] = (d.loc[:, ('pitzer', c)] + d.loc[:, ('MyAMI', c)]) / 2
    
    d.loc[:, ('csys_mid', 'H')] = 10**-d.loc[:, ('csys_mid', 'pHtot')]

    return d