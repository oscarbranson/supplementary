'''
Functions for using PHREEQC to calculate solution chemistry.

(c) 2019 : Oscar Branson : oscarbranson@gmail.com
'''
import os
import pkg_resources as pkgrs
import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod

def input_str_flex(param_dict,units='mol/kgw'):
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
    inp.append('units    {}'.format(units))
    for k, v in param_dict.items():
        if isinstance(v, str):
            inp.append('    {:20s}{:s}'.format(k, v))
        else:
            inp.append('    {:20s}{:.6e}'.format(k, v))
    
    post = """SELECTED_OUTPUT
    -pH                 true
    -temperature        true
    -alkalinity         true
    -ionic_strength     true
    -solution           false
    -sim                false
    -distance           false
    -time               false
    -step               false
    -state              false
    -totals Na Cl Ca C Li
    -m Cl- Na+ Ca+2 Li+
    -m OH- H+
    -m HCO3- CO3-2 CO2 H2CO3 # carbon
    -m NaCO3- NaHCO3 # sodium
    -m LiOH LiSO4
    -m CaHCO3+ CaOH+ CaCO3 # calcium
    -a Cl- Na+ Ca+2 Li+
    -a OH- H+
    -a HCO3- CO3-2 CO2 H2CO3
    -a NaCO3- NaHCO3 # sodium
    -a LiOH LiSO4
    -a CaHCO3+ CaOH+ CaCO3 # calcium
    -si Calcite Aragonite
END
    """
    
    return pre + '\n'.join(inp) + '\n' + post

def run_phreeqc(input_strings, dbase='pitzer', phreeq_path='/usr/local/lib/libiphreeqc.so'):
    """
    Run input string in phreeqc with specified database.

    Parameters
    ----------
    input_strings : str or list
        A valid phreeqc input string, or list of strings, with SELECTED_OUTPUT.
    dbase_path : str
        Path to valid phreeqc database (e.g. minteq.v4.dat)
    phreeq_path : str
        Path to iphreeqc shared library.

    Returns
    -------
    pandas.Series of calculated species
    """
    if not os.path.exists(dbase):
        dbase = get_local_dbase_path(dbase)
        if not os.path.exists(dbase):
            raise ValueError('PHREEQC database {} does not exist.'.format(dbase))
    
    if isinstance(input_strings, str):
        input_strings = {0: input_strings}

    phreeqc = phreeqc_mod.IPhreeqc(phreeq_path)
    phreeqc.load_database(dbase)
    outs = {}
    for k, input_string in input_strings.items():
        phreeqc.run_string(input_string)
        out = phreeqc.get_selected_output_array()
        outs[k] = dict(zip(out[0], out[1]))
    return outs

def get_local_dbase_path(dbase='pitzer'):
    """
    Returns the path to the phreeqc database packaged with this code.
    """
    return pkgrs.resource_filename('li_funks', 'db_phreeqc/' + dbase + '.dat')

def build_input_string(row):
    """
    Return a dict describing the solution from a dataframe row.
    """
    sol = dict(row.loc['solution_phreeq'])

    # _ = sol.pop('Li/Ca')
    # sol['pH'] = row.loc[('Conditions', 'pH')]
    # sol['temp'] = row.loc[('Conditions', 'TempC')]
    
    return input_str_flex(sol)
