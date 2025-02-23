def get_composition_after_calc(model, new_calcite_mol_cm2):
    m = model.copy()
    m.set_index('G_cum', inplace=True)
    return m.reindex(sorted([new_calcite_mol_cm2] + list(m.index))).interpolate().loc[new_calcite_mol_cm2, :]

def get_composition_after_hrs(model, hrs_to_calcify):
    m = model.copy()
    m.set_index('t_hr', inplace=True)
    return m.reindex(sorted([hrs_to_calcify] + list(m.index))).interpolate().loc[hrs_to_calcify, :]