'''
This script compares numerical results obtained from solving the radial Schrodinger equation,
(for charmonium and bottomonium systems) in Charm_and_Beauty.py
with experimental values from the Particle Data Group (PDG).
The computed and  experimental PDG values are displayed in side-by-side tables, 
including their absolute differences.

Results compared:
    - Energy levels
    - Wavefunction at origin: |Ψ(0)|²
    - Mean radius <r>

Tables are displayed using matplotlib, and data is structured with pandas.
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd  # Used for tabular data manipulation and formatting


def show_table_with_diff(computed, reference, headers, title, unit):
    '''
    Display a comparison table between computed and reference values.

    Parameters
    ----------
    computed : list of float
        Computed values (from numerical simulation).
    reference : list of float or None
        Experimental reference values (from PDG). Use None where unavailable.
    headers : list of str
        List containing the label for the computed column.
    title : str
        Title displayed above the table (e.g., "Charmonium Energies").
    unit : str
        Unit to be displayed (e.g., "MeV", "m", "m-³").

    Returns
    -------
    None. Displays a styled matplotlib table comparing values.
    '''
    rows = []

    for c, r in zip(computed, reference):
        # Handle missing experimental data
        if r is None:
            diff = "n/a"
            r_val = "n/a"
        else:
            diff = f"{abs(c - r):.4g} {unit}"  # Absolute difference
            r_val = r
        rows.append([c, r_val, diff])

    # Create DataFrame for table formatting and potential manipulation
    col_labels = [headers[0], f"Experimental PDG ({unit})", f"Difference ({unit})"]
    df = pd.DataFrame(rows, columns=col_labels).replace({np.nan: "n/a"})

    # Create and configure matplotlib figure
    fig, ax = plt.subplots(figsize=(9, len(df) * 0.6 + 0.8), dpi=150)
    ax.axis('off')  # Hide the axes frame

    # Add the data table
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        cellLoc='center',
        loc='center',
        colColours=["#dddddd"] * len(df.columns),  # Light gray header
    )
    table.scale(1.2, 1.3)  # Scale to improve legibility

    # Add title to the figure
    plt.subplots_adjust(top=0.85)
    plt.title(title, fontsize=13, weight='bold', pad=1)
    plt.show()


# ****** Experimental PDG values (source in report) ******
pdg_charmonium_s = [3097, 3686, 3773, 4039, 4191]
pdg_charmonium_p = [3415, 3511, 3556, 3525, None]
pdg_bottomonium_s = [9399, 9999, None, None, None]
pdg_bottomonium_p = [9859, 9893, 9912, 9899, None]

pdg_psi0_charmonium = [0.081, None, None, None, None]
pdg_psi0_bottomonium = [0.512, None, None, None, None]

pdg_r_charmonium = [0.25e-15, 0.45e-15, 0.65e-15, None, None]
pdg_r_bottomonium = [0.14e-15, 0.28e-15, 0.42e-15, None, None]


# ****** Computed values from Charm_and_Beauty.py ******

# All energies are in MeV, wavefunction density in m-³m-³, and radii in meters
computed_data = {
    "charmonium_s": {
        "energies": [983.37, 1381.36, 1829.10, 2359.75, 2940.15],
        "psi0": [2.83e48] * 5,
        "r_mean": [3.96e-51] * 5
    },
    "charmonium_p": {
        "energies": [979.22, 1368.92, 1824.95, 2343.17, 2936.01],
        "psi0": [0] * 5,  # â=1 â Ï(0) = 0 by angular momentum
        "r_mean": [3.96e-51] * 5
    },
    "bottomonium_s": {
        "energies": [3650.38, 5127.76, 6789.82, 8759.67, 10914.20],
        "psi0": [1.05e49] * 5,
        "r_mean": [1.07e-51] * 5
    },
    "bottomonium_p": {
        "energies": [3634.99, 5081.60, 6774.43, 8698.12, 10898.81],
        "psi0": [0] * 5,
        "r_mean": [1.07e-51] * 5
    }
}


# ****** Display Comparison Tables ******

# (1) Energies (S-wave and P-wave)
show_table_with_diff(computed_data["charmonium_s"]["energies"], pdg_charmonium_s,
                     ["Computed (MeV)"], "Charmonium Energies (l = 0)", "MeV")

show_table_with_diff(computed_data["charmonium_p"]["energies"], pdg_charmonium_p,
                     ["Computed (MeV)"], "Charmonium Energies (l = 1)", "MeV")

show_table_with_diff(computed_data["bottomonium_s"]["energies"], pdg_bottomonium_s,
                     ["Computed (MeV)"], "Bottomonium Energies (l = 0)", "MeV")

show_table_with_diff(computed_data["bottomonium_p"]["energies"], pdg_bottomonium_p,
                     ["Computed (MeV)"], "Bottomonium Energies (l = 1)", "MeV")

# (2) |Ψ(0)|² values
show_table_with_diff(computed_data["charmonium_s"]["psi0"], pdg_psi0_charmonium,
                     ["|Ψ(0)|² (m-³)"], "|Ψ(0)|² - Charmonium l = 0", "m-³")

show_table_with_diff(computed_data["charmonium_p"]["psi0"], [0]*5,
                     ["|Ψ(0)|² (m-³)"], "|Ψ(0)|² - Charmonium l = 1", "m-³")

show_table_with_diff(computed_data["bottomonium_s"]["psi0"], pdg_psi0_bottomonium,
                     ["|Ψ(0)|² (m-³)"], "|Ψ(0)|² - Bottomonium l = 0", "m-³")

show_table_with_diff(computed_data["bottomonium_p"]["psi0"], [0]*5,
                     ["|Ψ(0)|² (m-3)"], "|Ψ(0)|² - Bottomonium l = 1", "m-³")

# (3) Mean radius <r>
show_table_with_diff(computed_data["charmonium_s"]["r_mean"], pdg_r_charmonium,
                     ["<r>(m)"], "Mean Radius <r> - Charmonium l = 0", "m")

show_table_with_diff(computed_data["bottomonium_s"]["r_mean"], pdg_r_bottomonium,
                     ["<r>(m)"], "Mean Radius <r> - Bottomonium l = 0", "m")

