import pandas as pd
import numpy as np
import uncertainties.unumpy as unp


def raw_b_data(data_file, species=None):
    """
    Loads B/Ca data.
    
    Parameters
    ----------
    data_file : str
        Path to csv file containing data.
    species : str
        A string containing either "Genus species", "Genus" or "species",
        specifying which species to load.
"""
    d = pd.read_csv(data_file)
    ind = np.ones(d.shape[0], dtype=bool)
    
    if species is not None:
        sspl = species.split(' ')
        if len(sspl) == 2:
            genus = sspl[0]
            speci = sspl[1]
            ind = ind & (d.loc[:, 'Genus'] == genus) & (d.loc[:, 'Species'] == speci)
        elif len(sspl) == 1:
            genus = sspl[0]
            speci = sspl[0]
            ind = ind & ((d.loc[:, 'Genus'] == genus) | (d.loc[:, 'Species'] == speci))
        else:
            vspec = '\n    - ' + '\n    - '.join((d.loc[:, 'Genus'] + ' ' + d.loc[:, 'Species']).unique())
            raise ValueError('Dont understand species="{}". Please use a valid species name:'.format(species) + vspec)

    d.loc[:, 'who'] = [s[0] for s in d.Reference.str.split(',')]

    return d.loc[ind, :].apply(pd.to_numeric, errors='ignore')

def b_data(include=['This Study']):
    odat = pd.read_csv('data/B_compiled.csv', index_col=0, header=[0,1])

    odat.loc[:, ('Measured', 'KB')] = odat.loc[:, ('Measured', 'B/Caf')] / (odat.loc[:, ('Measured', '[B]sw')] / odat.loc[:, ('csys_mid', 'DIC')])

    # package uncertainties
    odat.loc[:, ('Measured', 'B/Caf')] = unp.uarray(odat.loc[:, ('Measured', 'B/Caf')], odat.loc[:, ('Uncertainties', 'se')])

    odat.loc[:, ('Measured', 'KB')] = odat.loc[:, ('Measured', 'B/Caf')] / (odat.loc[:, ('Measured', '[B]sw')] / odat.loc[:, ('csys_mid', 'DIC')])
    odat.loc[:, ('Measured', 'KB_HCO3')] = odat.loc[:, ('Measured', 'B/Caf')] / (odat.loc[:, ('Measured', '[B]sw')] / odat.loc[:, ('csys_mid', 'HCO3')])

    # a few derived parameters
    odat[('csys_mid', 'del_CO3')] = 1e6 * odat[('csys_mid', 'KspC')] / (1e-3 * odat[('Measured', '[Ca]sw')])
    odat[('csys_mid', 'BT_DIC')] = odat[('csys_mid', 'BT')] / odat[('csys_mid', 'DIC')]

    # comment out studies here to remove them
    # mdict = {
        # 'This Study': 'o', 
    #     'Haynes et al. (2017)': 's', 
    #     'Allen et al. (2011)': 'd',
#     'Howes et al. (2017)': '*'
# }

    if include is None:
        return odat
    else:
        return odat.loc[odat.loc[:, ('Measured', 'who')].isin(include)]