import numpy as np
import pandas as pd
import cbsyst as cb

from scipy.integrate import solve_ivp, cumulative_trapezoid
import warnings
warnings.simplefilter("ignore")


def constant_DIC(chamber_diameter_um, temperature, pHtot, J=None, pCO2=420, salinity=35, Ca_SW=10.2e-3, f_Ca=1, TE_SW=52e-3, F_M=0.0004, K_TE_G=2e-2, K_TE_J=0, G_k=1e-14, G_n=1.2, f_DIC=1, wall_thickness_um=3, hrs_to_run=12, method='Radau', rtol=1e-5, **solve_kwargs):
    """
    Simple transport model of calcifying fluid chemistry under constant DIC conditions.
    
    Parameters
    ----------
    chamber_diameter_um : float
        Diameter of the calcifying chamber in micrometers.
    temperature : float
        Temperature in degrees Celsius.
    pHtot : float
        Total scale pH of the seawater.
    J : float, optional
        Proton pumping rate in mol cm-2 s-1. If None, a temperature-dependent default is used.
    pCO2 : float, optional
        Partial pressure of CO2 in seawater in ppm. Default is 420 ppm. Used to calculate external water carbonate system.
    salinity : float, optional
        Salinity of the seawater. Default is 35. Used to calculate external water carbonate system.
    Ca_SW : float, optional
        Calcium concentration in seawater in mol kg-1. Default is 10.2ee-3 mol kg-1.
    f_Ca : float, optional
        Fraction of pumped protons that are balanced by calcium ions. Default is 1.
    TE_SW : float, optional
        Trace element concentration in seawater in mol kg-1. Default is 52e-3 mol kg-1 for Mg.
    F_M : float, optional
        Flushing rate of the calcifying fluid in s-1. Default is 0.0004 s-1 (40 minutes).
    K_TE_G : float, optional
        Partition coefficient for trace element during calcification. Default is 2e-2 for Mg.
    K_TE_J : float, optional
        Partition coefficient for trace element during Ca counter-transport associated with proton pumping. Default is 0.
    G_k : float, optional
        Pre-exponential factor for calcification rate. Default is 1e-14 mol cm-2 s-1.
    G_n : float, optional
        Exponent for saturation state dependence of calcification rate. Default is 1.2.
    f_DIC : float, optional
        Factor by which DIC is elevated in the calcifying fluid relative to seawater. Default is 1.
    wall_thickness_um : float, optional
        Thickness of the calcifying chamber wall in micrometers. Default is 3 um.
    hrs_to_run : float, optional
        Duration of the simulation in hours. Default is 12 hours.
    method : str, optional
        Integration method for solve_ivp. Default is 'Radau'.
    rtol : float, optional
        Relative tolerance for the ODE solver. Default is 1e-5.
    **solve_kwargs : additional keyword arguments
        Additional arguments passed to solve_ivp.
    """
    
    initial_csys = cb.Csys(pHtot=pHtot, pCO2=pCO2, S_in=salinity, T_in=temperature, unit='mol')
    DIC = initial_csys.DIC.item() * f_DIC

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


    sol = solve_ivp(transport_model, (0, hrs_to_run * 60 * 60), [initial_csys.TA, initial_csys.Ca, TE_SW], method=method, rtol=rtol, **solve_kwargs)
    
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


def constant_DIC_rate(chamber_diameter_um, temperature, pHtot, J=None, pCO2=420, salinity=35, Ca_SW=10.2e-3, f_Ca=1,TE_SW=52e-3, F_M=0.0004, K_TE_G=2e-2, K_TE_J=0, G_k=1e-14, G_n=1.2, G_mid=2e-9/2, f_DIC=1, wall_thickness_um=3, hrs_to_run=12, method='Radau', rtol=1e-5, **solve_kwargs):
    
    initial_csys = cb.Csys(pHtot=pHtot, pCO2=pCO2, S_in=salinity, T_in=temperature, unit='mol')
    DIC = initial_csys.DIC.item() * f_DIC

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


    sol = solve_ivp(transport_model, (0, hrs_to_run * 60 * 60), [initial_csys.TA, initial_csys.Ca, TE_SW], method=method, rtol=rtol, **solve_kwargs)
    
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
