import numpy as np
import matplotlib.pyplot as plt

def time_resolved(output, vars=None, hrs_to_calcify=None, new_calcite_mol_cm2=None, axs=None, label=None, **kwargs):
    
    if vars is None:
        vars = [c for c in output.columns if c not in ['t', 't_hr']]
        
    if axs is None:
        fig, axs = plt.subplots(len(vars), 1, figsize=(5, 1.3*len(vars)), sharex=True, constrained_layout=True)
    else:
        fig = axs[0].get_figure()
        
    for ax, var in zip(axs, vars):
        if 'TE/Ca' in var:
            y = np.log(output[var])
            ax.set_ylabel('log(' + var + ')')
        else:
            y = output[var]
            ax.set_ylabel(var)
        
        ax.plot(output.t / 3600, y, label=label, **kwargs)
        
        if hrs_to_calcify is not None:
            ax.scatter(hrs_to_calcify, np.interp(hrs_to_calcify, output.t_hr, y), color='k')
            # ax.axvline(hrs_to_calcify)
                
        # if var == 'G_cum':
            # if new_calcite_mol_cm2 is not None:
                # ax.scatter(np.interp(new_calcite_mol_cm2, output.G_cum, output.t_hr), new_calcite_mol_cm2, color='k')
                # ax.axhline(new_calcite_mol_cm2)
            
    return fig, axs