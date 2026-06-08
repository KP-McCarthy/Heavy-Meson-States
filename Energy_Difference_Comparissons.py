'''
This script compares computed energy levels with experimental data from the PDG.
For each of the four systems (charmonium and bottomonium, S- and P-waves), it:

    (1) Plots the computed vs PDG energy levels.
    (2) Calculates and visualizes differences in level spacings (residuals).
    
Each comparison is shown in a two-part figure:
    - Top subplot: energy levels
    - Bottom subplot: spacing residuals (ΔE_computed - ΔE_PDG)
    
Missing data points (None) are automatically filtered.
'''

import matplotlib.pyplot as plt
import numpy as np


# ****** FUNCTION: PLOT ENERGY COMPARISON ******

def plot_energy_comparison(n_values, computed, pdg, label):
    '''
    Plots a 2-panel comparison between computed and PDG energy levels.

    Parameters
    ----------
    n_values : list of int
        Indices for the energy levels (e.g., 0 to 4).
    computed : list of float
        Computed energy levels (MeV).
    pdg : list of float
        Experimental PDG energy levels (MeV).
    label : str
        Label describing the system (e.g., "Charmonium ℓ = 0").

    Returns
    -------
    None. Displays the plot.
    '''
    
    # Set up a 2-row figure: top = energies, bottom = spacing residuals
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                            gridspec_kw={'height_ratios': [3, 1]})
    
    # ****** PLOT 1: ENERGY LEVELS ******
    axs[0].plot(n_values, computed, 'o-', label='Computed', color='blue')
    axs[0].plot(n_values, pdg, 'o--', label='PDG', color='orange')
    axs[0].set_ylabel("Energy (MeV)")
    axs[0].set_title(f"{label} — Energy Levels")
    axs[0].legend()
    axs[0].grid(True)

    # ****** PLOT 2: RESIDUALS OF ENERGY SPACING ******
    spacing_computed = np.diff(computed)
    spacing_pdg = np.diff(pdg)

    # Residuals: difference in spacings (ΔE_computed - ΔE_PDG)
    spacing_diff = spacing_computed - spacing_pdg
    n_spacing = [f"{i}-{i+1}" for i in range(len(spacing_diff))]

    axs[1].bar(n_spacing, spacing_diff, color='purple')
    axs[1].axhline(0, color='black', linestyle='--')
    axs[1].set_ylabel("ΔΔE (MeV)")
    axs[1].set_xlabel("State spacing")
    axs[1].set_title("Residuals: Computed ΔE - PDG ΔE")
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()



# ****** SYSTEM DATA ******

# Experimental and computed energy levels (in MeV)

# Charmonium ℓ = 0 (S-wave)
pdg_charmonium_s = [3097, 3686, 3773, 4039, 4191]
computed_charmonium_s = [983.37, 1381.36, 1829.10, 2359.75, 2940.15]

# Charmonium ℓ = 1 (P-wave)
pdg_charmonium_p = [3415, 3511, 3556, 3525, None]  # Final level missing
computed_charmonium_p = [979.22, 1368.92, 1824.95, 2343.17, 2936.01]

# Bottomonium ℓ = 0 (S-wave)
pdg_bottomonium_s = [9399, 9999, 10355, 10579, 10860]
computed_bottomonium_s = [3650.38, 5127.76, 6789.82, 8759.67, 10914.20]

# Bottomonium ℓ = 1 (P-wave)
pdg_bottomonium_p = [9859, 9893, 9912, 9899, None]  # Final level missing
computed_bottomonium_p = [3634.99, 5081.60, 6774.43, 8698.12, 10898.81]



# ****** DATA CLEANING FUNCTION ******

def clean_data(computed, pdg):
    '''
    Removes any levels where either computed or PDG value is missing (None).

    Parameters
    ----------
    computed : list
        Computed energy values.
    pdg : list
        Experimental PDG values.

    Returns
    -------
    tuple:
        cleaned_computed : list of valid computed values
        cleaned_pdg : list of corresponding PDG values
        indices : list of int indices (for x-axis plotting)
    '''
    cleaned_c, cleaned_p = [], []
    for c, p in zip(computed, pdg):
        if c is not None and p is not None:
            cleaned_c.append(c)
            cleaned_p.append(p)
    return cleaned_c, cleaned_p, list(range(len(cleaned_c)))



# ****** PLOT ALL FOUR SYSTEMS ******

# Charmonium ℓ = 0
c_c, p_c, n_vals = clean_data(computed_charmonium_s, pdg_charmonium_s)
plot_energy_comparison(n_vals, c_c, p_c, "Charmonium ℓ = 0")

# Charmonium ℓ = 1
c_c, p_c, n_vals = clean_data(computed_charmonium_p, pdg_charmonium_p)
plot_energy_comparison(n_vals, c_c, p_c, "Charmonium ℓ = 1")

# Bottomonium ℓ = 0
c_c, p_c, n_vals = clean_data(computed_bottomonium_s, pdg_bottomonium_s)
plot_energy_comparison(n_vals, c_c, p_c, "Bottomonium ℓ = 0")

# Bottomonium ℓ = 1
c_c, p_c, n_vals = clean_data(computed_bottomonium_p, pdg_bottomonium_p)
plot_energy_comparison(n_vals, c_c, p_c, "Bottomonium ℓ = 1")
