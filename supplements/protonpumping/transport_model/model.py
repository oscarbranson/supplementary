import numpy as np
import pandas as pd
import cbsyst as cb

from scipy.integrate import solve_ivp, cumulative_trapezoid
import warnings
warnings.simplefilter("ignore")


def constant_DIC(chamber_diameter_um, temperature, pHtot, J=None, pCO2=420, salinity=35, Ca_SW=10.2e-3, f_Ca=1, TE_SW=52e-3, F_M=0.0004, K_TE_G=2e-2, K_TE_J=0, G_k=1e-14, G_n=1.2, wall_thickness_um=3, hrs_to_run=12, method='Radau', rtol=1e-5, **solve_kwargs):
    
    initial_csys = cb.Csys(pHtot=pHtot, pCO2=pCO2, S_in=salinity, T_in=temperature, unit='mol')
    DIC = initial_csys.DIC.item()

    chamber_radius_cm = chamber_diameter_um / 2 * 1e-4  # cm
    chamber_surface_area_cm2 = 4 * np.pi * chamber_radius_cm**2 / 2  # cm2
    
    TEST_THICKNESS_CM = wall_thickness_um * 1e-4
    TEST_POROSITY = 0.25
    
    calc_vol_L = (1 - TEST_POROSITY) * 2. / 3 * np.pi * ((chamber_radius_cm + TEST_THICKNESS_CM/2)**3 - (chamber_radius_cm - TEST_THICKNESS_CM/2)**3) * 1e-3
    
    # calc_vol_L = (chamber_surface_area_cm2 * 0.0003) / 1e3  # L, assuming 3um test thickness

    CF_mass_total = calc_vol_L * 1.025  # mass of calcifying fluid, kg
    M_CF = CF_mass_total / chamber_surface_area_cm2  # kg cm-2

    calc_mass_kg = calc_vol_L * 2.71  # total mass of calcite, kg
    calc_moles = calc_mass_kg * 1e3 / 100.09  # moles of calcite
    C_CF = calc_moles / chamber_surface_area_cm2  # moles cm-2
    
    if J is None:
        J = np.polyval([0.24378071, -2.4807883], temperature) * 1e-9  # mol cm-2 s-1, parameters from data
    
    def transport_model(t, state):
        TA, Ca, TE = state
        
        # csys = cb.Csys(TA=TA, DIC=initial_csys.DIC, S_in=salinity, T_in=temperature, unit='mol')
        # OmegaC = csys.CO3.item() * Ca / csys.Ks.KspC.item()

        # approximate carbon system as CO3 = TA - DIC
        CO3 = min(TA - DIC, DIC)
        OmegaC = CO3 * Ca / initial_csys.Ks.KspC.item()
        # print(csys.CO3.item(), CO3, TA, initial_csys.DIC.item())
        # print(csys.CO3.item(), Ca, csys.Ks.KspC.item(), OmegaC)
            
        G = G_k * max(0, (OmegaC - 1))**G_n
                
        dCa = (
            0.5 * J * f_Ca / M_CF -  # pumping
            G / M_CF +  # calcification
            F_M * (Ca_SW - Ca)  # seawater_exchange
            )
        dTA = (
            J / M_CF -  # pumping
            2 * G / M_CF +  # calcification
            F_M * (initial_csys.TA - TA))  # seawater_exchange
        dTE = (
            F_M * (TE_SW - TE) -  # seawater_exchange
            K_TE_G * TE / Ca * G / M_CF +  # calcification
            K_TE_J * TE_SW / Ca_SW * J / M_CF  # pumping
            )
        
        return [dTA.item(), dCa, dTE]


    sol = solve_ivp(transport_model, (0, hrs_to_run * 60 * 60), [initial_csys.TA.item(), initial_csys.Ca.item(), TE_SW], method=method, rtol=rtol, **solve_kwargs)
    
    out = {k: v for k, v in zip(['TA', 'Ca', 'TE'], sol.y)}
    out['t'] = sol.t
    out['t_hr'] = sol.t / 3600
    
    out_csys = cb.Csys(TA=out['TA'], DIC=initial_csys.DIC, S_in=salinity, T_in=temperature, unit='mol')

    out['OmegaC'] = out_csys.CO3 * out['Ca'] / out_csys.Ks.KspC.item()
    out['G'] = G_k * (out['OmegaC'] - 1)**G_n  # mol cm-2 s-1
    
    out['G_cum'] = np.concatenate([[0], cumulative_trapezoid(out['G'], out['t'])])  # mol cm-2
    
    out = pd.DataFrame.from_dict(out)
    
    # calculate test composition
    out['fluid_TE/Ca'] = out.TE / out.Ca
    out['mineral_TE/Ca'] = K_TE_G * out['fluid_TE/Ca']

    out['test_TE/Ca'] = (out['mineral_TE/Ca'] * out.G_cum.diff()).cumsum() / out.G_cum
    
    out['temperature'] = temperature
    
    return out


def constant_DIC_rate(chamber_diameter_um, temperature, pHtot, J=None, pCO2=420, salinity=35, Ca_SW=10.2e-3, f_Ca=1,TE_SW=52e-3, F_M=0.0004, K_TE_G=2e-2, K_TE_J=0, G_k=1e-14, G_n=1.2, G_mid=2e-9/2, wall_thickness_um=3, hrs_to_run=12, method='Radau', rtol=1e-5, **solve_kwargs):
    
    initial_csys = cb.Csys(pHtot=pHtot, pCO2=pCO2, S_in=salinity, T_in=temperature, unit='mol')
    DIC = initial_csys.DIC.item()

    chamber_radius_cm = chamber_diameter_um / 2 * 1e-4  # cm
    chamber_surface_area_cm2 = 4 * np.pi * chamber_radius_cm**2 / 2  # cm2
    
    TEST_THICKNESS_CM = wall_thickness_um * 1e-4
    TEST_POROSITY = 0.25
    
    calc_vol_L = (1 - TEST_POROSITY) * 2. / 3 * np.pi * ((chamber_radius_cm + TEST_THICKNESS_CM/2)**3 - (chamber_radius_cm - TEST_THICKNESS_CM/2)**3) * 1e-3
    
    # calc_vol_L = (chamber_surface_area_cm2 * 0.0003) / 1e3  # L, assuming 3um test thickness

    CF_mass_total = calc_vol_L * 1.025  # mass of calcifying fluid, kg
    M_CF = CF_mass_total / chamber_surface_area_cm2  # kg cm-2

    calc_mass_kg = calc_vol_L * 2.71  # total mass of calcite, kg
    calc_moles = calc_mass_kg * 1e3 / 100.09  # moles of calcite

    C_CF = calc_moles / chamber_surface_area_cm2  # moles cm-2
    
    if J is None:
        J = np.polyval([0.24378071, -2.4807883], temperature) * 1e-9  # mol cm-2 s-1, parameters from data
    
    def transport_model(t, state):
        TA, Ca, TE = state
        
        # csys = cb.Csys(TA=TA, DIC=initial_csys.DIC, S_in=salinity, T_in=temperature, unit='mol')
        # OmegaC = csys.CO3.item() * Ca / csys.Ks.KspC.item()
        # print(csys.CO3.item(), Ca, csys.Ks.KspC.item(), OmegaC)
        
        # approximate carbon system as CO3 = TA - DIC
        CO3 = min(TA - DIC, DIC)
        OmegaC = CO3 * Ca / initial_csys.Ks.KspC.item()

        G = G_k * max(0, (OmegaC - 1))**G_n
        
        # print(state, csys.OmegaC, G, J)
        
        dCa = (
            0.5 * J * f_Ca / M_CF -  # pumping
            G / M_CF +  # calcification
            F_M * (Ca_SW - Ca)  # seawater_exchange
            )
        dTA = (
            J / M_CF -  # pumping
            2 * G / M_CF +  # calcification
            F_M * (initial_csys.TA - TA))  # seawater_exchange
        dTE = (
            F_M * (TE_SW - TE) -  # seawater_exchange
            (K_TE_G * G / G_mid) * TE / Ca * G / M_CF +  # calcification
            K_TE_J * TE_SW / Ca_SW * J / M_CF  # pumping
            )
        
        return [dTA.item(), dCa, dTE]


    sol = solve_ivp(transport_model, (0, hrs_to_run * 60 * 60), [initial_csys.TA.item(), initial_csys.Ca.item(), TE_SW], method=method, rtol=rtol, **solve_kwargs)
    
    out = {k: v for k, v in zip(['TA', 'Ca', 'TE'], sol.y)}
    out['t'] = sol.t
    out['t_hr'] = sol.t / 3600
    
    out_csys = cb.Csys(TA=out['TA'], DIC=initial_csys.DIC, S_in=salinity, T_in=temperature, unit='mol')

    out['OmegaC'] = out_csys.CO3 * out['Ca'] / out_csys.Ks.KspC.item()
    out['G'] = G_k * (out['OmegaC'] - 1)**G_n
    
    out['G_cum'] = np.concatenate([[0], cumulative_trapezoid(out['G'], out['t'])])
    
    out = pd.DataFrame.from_dict(out)
    
    # calculate test composition
    out['fluid_TE/Ca'] = out.TE / out.Ca
    out['mineral_TE/Ca'] = K_TE_G * out['fluid_TE/Ca']

    out['test_TE/Ca'] = (out['mineral_TE/Ca'] * out.G_cum.diff()).cumsum() / out.G_cum
    
    out['temperature'] = temperature
    
    return out
