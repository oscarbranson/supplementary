"""
Functions for calculating C speciation in seawater using PHREEQC
"""

import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod

def input_str_flex(param_dict):
    """
    Generate PHREEQC solution input string from dict

    Parameters
    ----------
    param_dict : dict
        key: value pairs that are valid entries in a
        PHREEQC solution definition. (See
        https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-56.html)
    
    Returns
    -------
    str : input string for PHREEQC
    """
    pre = "SOLUTION 1\n"
    
    inp = []
    for k, v in param_dict.items():
        if isinstance(v, str):
            inp.append('    {:20s}{:s}'.format(k, v))
        else:
            inp.append('    {:20s}{:.8f}'.format(k, v))
    
    post = """SELECTED_OUTPUT
    -pH
    -temperature
    -alkalinity
    -ionic_strength
    -totals Cl Na Mg K B Ca C S(6) Br
    -m OH- H+
    # pitzer outputs
    -m B(OH)4- B(OH)3 CaB(OH)4+ B3O3(OH)4- B4O5(OH)4-2  # boron
    -m HCO3- CO3-2 CO2 H2CO3 # carbon
    -m HSO4- SO4-2  # sulphate
    -si Calcite Aragonite
END
    """
    
    return pre + '\n'.join(inp) + '\n' + post

def run_phreeqc(input_string, dbase_path, phreeq_path='/usr/local/lib/libiphreeqc.so'):
    """
    Run input string in phreeqc with specified database.

    Parameters
    ----------
    input_string : str
        Valid phreeqc input string with SELECTED_OUTPUT.
    dbase_path : str
        Path to valid phreeqc database (e.g. minteq.v4.dat)
    phreeq_path : str
        Path to iphreeqc shared library.

    Returns
    -------
    pandas.Series of calculated species
    """
    phreeqc = phreeqc_mod.IPhreeqc(phreeq_path)
    phreeqc.load_database(dbase_path)
    phreeqc.run_string(input_string)
    out = phreeqc.get_selected_output_array()
    return dict(zip(out[0], out[1]))

def Kcalc(reactant_A, reactant_B, product):
    """
    Calculate equilibrium constant from reactant activities.

    Parameters
    ----------
    reactant_A, reactant_B : float or array-like
        Activities of reactants in solution
    product : float or array-like
        Activity of product in solution.

    Returns
    -------
    float or array-like : Equilibrium Constant (K)
    """
    return reactant_A * reactant_B / product

def Ks_from_pitzer(pitzer_output):
    Ks = {
        'K1': Kcalc(pitzer_output['m_H+(mol/kgw)'], pitzer_output['m_HCO3-(mol/kgw)'], 
                    pitzer_output['m_H2CO3(mol/kgw)'] + pitzer_output['m_CO2(mol/kgw)']),
        'K2': Kcalc(pitzer_output['m_H+(mol/kgw)'], pitzer_output['m_CO3-2(mol/kgw)'], pitzer_output['m_HCO3-(mol/kgw)']),
        'KB': Kcalc(pitzer_output['m_H+(mol/kgw)'], pitzer_output['m_B(OH)4-(mol/kgw)'], pitzer_output['m_B(OH)3(mol/kgw)']),
        'KS': Kcalc(pitzer_output['m_H+(mol/kgw)'], pitzer_output['m_SO4-2(mol/kgw)'], pitzer_output['m_HSO4-(mol/kgw)']),
        'KW': Kcalc(pitzer_output['m_H+(mol/kgw)'], pitzer_output['m_OH-(mol/kgw)'], 1),
        }
    return Ks

def calc_pitzer_Ks(sw, database="./mg_funks/db_phreeqc/pitzer.dat", phreeq_path='/usr/local/lib/libiphreeqc.so'):
    inp = input_str_flex(sw)
    c = run_phreeqc(inp, database, phreeq_path)
    
    return Ks_from_pitzer(c)

def calc_pitzer_fKs(sw, cond, database="./mg_funks/db_phreeqc/pitzer.dat", phreeq_path='/usr/local/lib/libiphreeqc.so'):
    swinp = input_str_flex(sw)
    swc = run_phreeqc(swinp, database, phreeq_path)
    swKs = Ks_from_pitzer(swc)

    condinp = input_str_flex(cond)
    condc = run_phreeqc(condinp, database, phreeq_path)
    condKs = Ks_from_pitzer(condc)

    fKs = {}
    for k in swKs.keys():
        fKs[k] = condKs[k] / swKs[k]
    
    return fKs
