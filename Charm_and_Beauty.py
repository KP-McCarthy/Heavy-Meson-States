'''

The following script is for project title:
    "The Potential of charm and Beauty"
By Kevin McCarthy and Daniel Farrell

Solving the radial Schrodinger Equation using a Cornell potential for:
    (1) Charmonium: S-waves and P-waves
    (2) Bottomonium: S-waves and P-waves

For each of the 4 systems, for lowest 5 eigenenergies:
- Plotting graphs of:
    (1) Bracket finding
    (2) Energy Spectrum
    (3) Normalsied Wavefunction versus rho and r
    
- Outputting to console:
    (1) Eigenvalues
    (2) Wavefunction squared at zero
    (3) Mean radius


Glossary of Parts;
(1) PART 1: Constants and 4 Systems Set Up
    ****** CONTROL LOOP FOR SYSTEMS ******
(2) PART 2: Potential and ODE system
(3) PART 3: Shooting Method
(4) PART 4: Find Brackets and Plot Process
(5) PART 5: Bisection Method
(6) PART 6: Find Eigenvalues and Plot (Dimensionless)
(7) PART 7: Solve and Cache All Wavefunctions Once
(8) PART 8: Solve Schrodinger Equation, Normalise, and Plot (Dimensionless)
(9) PART 9: |Ψ(0)|² and mean radius
(10) PART 10: Conversion to Standard Units
    ****** END OF CONTROL LOOP FOR SYSTEMS ******
(11) PART 11: Wavefunction Intensity Heatmap

End of Intro
'''


# Importing Libraries 

import numpy as np
import matplotlib.pyplot as plt

# ****** PART 1: Constants and 4 Systems Set Up ******

# alpha prime is written as just alpha for brevity, tidyness
alpha = 0.472 # Not actually dimensionless, see dimensionless potential function
h = 0.001 # step size for rk4 and rho step in shooting function

rho0 = 1e-5 # initial close to zero rho
rhomax = 20 # chosen to include 5 lowest eigenenrgies

hbar_SI = 1.05e-34       # SI units (J·s), used in potential
hbar_ev = 6.5821e-16     # eV·s, used for conversions to meters

c = 299792458  # m/s

# Define system combinations
systems = [
    ('Charmonium', 1.32e9, 0),  # S-wave
    ('Charmonium', 1.32e9, 1),  # P-wave
    ('Bottomonium', 4.9e9, 0),  # S-wave
    ('Bottomonium', 4.9e9, 1)   # P-wave
]

# Dictionary to hold data for all 4 systems
all_wavefunction_data = {}


# ****** END OF PART 1: Constants and 4 Systems Set Up ******


# ****** CONTROL LOOP FOR SYSTEMS ******

# Iterates all solving for each of the 4 systems
for system_name, m, ell in systems:
    print(f"\n\n Solving for {system_name}, ℓ = {ell} \n")

    # Derived quantities per system
    m2c4 = (m**2)*(c**4)  # For use in kappa, avoids runtime error 
    kappa = (440**2) / m2c4

    # Reset wavefunction cache per system
    wavefunction_data = {}
    
    # Defining function for initial conditions
    def get_initial_conditions(ell):
        '''
        Returns initial conditions based on angular momentum.
        For S-wave (ell=0): small ψ(0), zero derivative.
        For P-wave (ell=1): zero ψ(0), small derivative.
        '''
        return (1e-8, 0) if ell == 0 else (0, 1e-8)

# ****** ALL FOLLOWING CODE IS WITHIN CONTROL LOOP FOR SYSTEMS ******


    # ****** PART 2: Potential and ODE system ******
    
    # Defining the Cornell Potential in dimensionless form
    def V_rho(rho, alpha, kappa):
        '''
        Compute the dimensionless Cornell potential at a given rho.
    
        Parameters
        ----------
        rho : float
            Dimensionless radial coordinate.
        alpha : float
            Effective coupling constant (dimensionless).
        kappa : float
            Linear confinement coefficient (dimensionless).
    
        Returns
        -------
        float
            Value of the potential V(rho).
        '''
        if rho == 0:
            return 0  # avoid singularity
        return -(alpha / rho * (hbar_SI**2) * (c**2)) + (kappa * rho)
    
    
    # Defining Schrodinger equation in Dimensionless form
    def schrodinger(rho, y, params):
        '''
        Evaluate the right-hand side of the radial Schrödinger equation.
    
        Parameters
        ----------
        rho : float
            Dimensionless radial coordinate.
        y : ndarray
            Array containing [ψ, ψ'] at current rho.
        params : tuple
            Tuple containing (epsilon, alpha, kappa, ell) for the potential and angular momentum.
    
        Returns
        -------
        ndarray
            Array of derivatives [dψ/drho, d²ψ/drho²].
        '''
        psi, phi = y
        epsilon, alpha, kappa, ell = params
        if rho == 0:
            rho = 1e-8  # avoid division by zero
        V = V_rho(rho, alpha, kappa) # tidyness, for use in schrondinger calculation
        dpsi_drho = phi
        dphi_drho = (ell*(ell+1)/rho**2 + V - epsilon) * psi
        return np.array([dpsi_drho, dphi_drho])
    
    
    # 4th Order Runge-Katta
    def rk4_step(f, rho, y, h, params):
        '''
        Perform a single 4th-order Runge-Kutta integration step.
    
        Parameters
        ----------
        f : function
            The function defining the system of ODEs, must return dy/dx.
        rho : float
            The current value of the independent variable.
        y : ndarray
            Current values of the dependent variables.
        h : float
            Step size for integration.
        params : tuple
            Parameters to pass to the ODE function.
    
        Returns
        -------
        ndarray
            Updated values of y after one RK4 step.
        '''
        k1 = f(rho, y, params)
        k2 = f(rho + 0.5*h, y + 0.5*h*k1, params)
        k3 = f(rho + 0.5*h, y + 0.5*h*k2, params)
        k4 = f(rho + h, y + h*k3, params)
        return y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    
    # ****** END OF PART 2: Potential and ODE system ******
    
    # ****** PART 3: Shooting Method ******
    
    def shoot(epsilon, alpha, kappa, ell, rho0, rhomax, h):
        '''
        Integrate the Schrödinger equation from rho0 to rhomax using the shooting method.
    
        Parameters
        ----------
        epsilon : float
            Trial energy eigenvalue.
        alpha : float
            Coupling constant for the potential.
        kappa : float
            Confinement term coefficient.
        ell : int
            Angular momentum quantum number.
        rho0 : float
            Starting value for rho (close to zero).
        rhomax : float
            Upper limit of integration.
        h : float
            Step size.
    
        Returns
        -------
        float
            Value of the wavefunction ψ at rho = rhomax.
        '''
        params = (epsilon, alpha, kappa, ell)
        # setting initial conditions 
        psi0, dpsi0 = get_initial_conditions(ell)

        y = np.array([psi0, dpsi0])
    
        rho = rho0
        while rho < rhomax:
            y = rk4_step(schrodinger, rho, y, h, params)
            rho += h
    
        return y[0]  # return ψ at rhomax
    
    
    # ****** END OF PART 3: Shooting Method ******
    
    # ****** PART 4: Find Brackets and Plot Process ******
    
    def find_brackets_plot(alpha, kappa, rho0, rhomax, h, eps_min, eps_max, eps_steps):
        '''
        Identify intervals of epsilon with sign changes in ψ(rho_max), indicating roots.
    
        Parameters
        ----------
        alpha : float
            Coupling constant for the potential.
        kappa : float
            Confinement term coefficient.
        rho0 : float
            Initial radial value.
        rhomax : float
            Final radial value.
        h : float
            Step size.
        eps_min : float
            Minimum trial epsilon.
        eps_max : float
            Maximum trial epsilon.
        eps_steps : int
            Number of trial steps between eps_min and eps_max.
    
        Returns
        -------
        list of tuples
            Each tuple contains (lower_bound, upper_bound) bracketing an eigenvalue.
        '''
        epsilons = np.linspace(eps_min, eps_max, eps_steps)
        values = []
        brackets = []
    
        prev_val = shoot(epsilons[0], alpha, kappa, ell, rho0, rhomax, h)
        values.append(prev_val)
    
        for eps in epsilons[1:]:
            val = shoot(eps, alpha, kappa, ell, rho0, rhomax, h)
            values.append(val)
            if prev_val * val < 0:
                brackets.append((eps - (epsilons[1] - epsilons[0]), eps))
            prev_val = val
    
        plt.plot(epsilons, values)
        plt.axhline(0, color='black', linestyle='--')
        plt.xlabel('ε')
        plt.ylabel('ψ(ρ_max)')
        plt.title('Finding sign changes: ψ(ρ_max) vs ε')
        plt.grid(True)
        plt.show()
    
        return brackets
    
    
    # ****** END OF PART 4: Find Brackets and Plot Process ******
    
    # ****** PART 5: Bisection Method ******
    
    def bisection(alpha, kappa, ell, rho0, rhomax, h, a, b, tol=1e-6, max_iter=100):
        '''
        Find an eigenvalue using the bisection method within a bracket.
    
        Parameters
        ----------
        alpha : float
            Coupling constant for the potential.
        kappa : float
            Confinement term coefficient.
        ell : int
            Angular momentum quantum number.
        rho0 : float
            Starting value of radial coordinate.
        rhomax : float
            End value of radial coordinate.
        h : float
            Step size for integration.
        a : float
            Lower bound of epsilon bracket.
        b : float
            Upper bound of epsilon bracket.
        tol : float, optional
            Tolerance for stopping, by default 1e-6.
        max_iter : int, optional
            Maximum number of iterations, by default 100.
    
        Returns
        -------
        float
            Estimated eigenvalue ε that satisfies boundary conditions.
        '''
        fa = shoot(a, alpha, kappa, ell, rho0, rhomax, h)
        fb = shoot(b, alpha, kappa, ell, rho0, rhomax, h)
    
        if fa * fb > 0:
            raise ValueError("No sign change in bracket!")
    
        for _ in range(max_iter):
            c = (a + b) / 2
            fc = shoot(c, alpha, kappa, ell, rho0, rhomax, h)
    
            if abs(fc) < tol:
                return c
    
            if fa * fc < 0:
                b, fb = c, fc
            else:
                a, fa = c, fc
    
        return (a + b) / 2  # best guess
    
    
    # ****** END OF PART 5: Bisection Method ******
    
    # ****** PART 6: Find Eigenvalues and Plot ******
    
    brackets = find_brackets_plot(alpha, kappa, rho0, rhomax, h, 0.5, 3.0, 200)
    
    eigenvalues = []
    for a, b in brackets[:5]:  # first 5 eigenstates
        eigen = bisection(alpha, kappa, ell, rho0, rhomax, h, a, b)
        eigenvalues.append(eigen)
        print(f"Eigenvalue found: ε = {eigen:.6f}")
    
    plt.figure()
    for n, ev in enumerate(eigenvalues):
        En = ev*m # mass is already in eV, no need to multiply by c^2
        plt.hlines(En, xmin=0, xmax=1, label=f'n={n}')
    plt.title(f'Energy spectrum: {system_name}, ℓ = {ell}')
    plt.xlabel('State index (arbitrary)')
    plt.ylabel('ε (dimensionless energy)')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # ****** END OF PART 6: Find Eigenvalues and Plot ******
    
    # ****** PART 7: Solve and Cache All Wavefunctions Once ******
    
    # This block prevents repeated calls to solve_wavefunction and normalize_wavefunction.
   
    def solve_wavefunction(epsilon, alpha, kappa, ell, rho0, rhomax, h):
        '''
        Integrate the radial Schrödinger equation to obtain the wavefunction.
    
        Parameters
        ----------
        epsilon : float
            Energy eigenvalue.
        alpha : float
            Coupling constant.
        kappa : float
            Confinement strength.
        ell : int
            Angular momentum quantum number.
        rho0 : float
            Initial rho value.
        rhomax : float
            Maximum rho to integrate to.
        h : float
            Step size.
    
        Returns
        -------
        tuple of ndarrays
            rho_values : Array of rho values.
            psi_values : Array of wavefunction values at corresponding rho.
        '''
        params = (epsilon, alpha, kappa, ell)
        
        # setting initial conditions 
        psi0, dpsi0 = get_initial_conditions(ell)
        
        y = np.array([psi0, dpsi0])
    
        rho_values = [rho0]
        psi_values = [psi0]
    
        rho = rho0
        while rho < rhomax:
            y = rk4_step(schrodinger, rho, y, h, params)
            rho += h
            rho_values.append(rho)
            psi_values.append(y[0])
    
        return np.array(rho_values), np.array(psi_values)
    
    
    def normalize_wavefunction(rho_values, psi_values):
        '''
        Normalize a wavefunction using the spherical volume element.
    
        Parameters
        ----------
        rho_values : ndarray
            Array of dimensionless radial coordinates.
        psi_values : ndarray
            Unnormalized wavefunction values.
    
        Returns
        -------
        ndarray
            Normalized wavefunction.
        '''
        n = len(rho_values)
        integral = 0
        for i in range(n-1):
            drho = rho_values[i+1] - rho_values[i]
            integrand_i = rho_values[i]**2 * psi_values[i]**2
            integrand_ip1 = rho_values[i+1]**2 * psi_values[i+1]**2
            integral += 0.5 * (integrand_i + integrand_ip1) * drho
    
        norm_factor = np.sqrt(integral)
        normalized_psi = psi_values / norm_factor
        return normalized_psi
    
    # Precompute and store wavefunctions to avoid recomputation
    wavefunction_data = {}
    
    for i, eigenvalue in enumerate(eigenvalues):
        print(f"Caching wavefunction for state n={i} with ε={eigenvalue:.6f}...")
        rho_vals, psi_vals = solve_wavefunction(eigenvalue, alpha, kappa, ell, rho0, rhomax, h)
        norm_psi = normalize_wavefunction(rho_vals, psi_vals)
        r_vals = rho_vals * (hbar_ev / (m * c))  # Convert to r (meters)
    
        wavefunction_data[i] = {
            "epsilon": eigenvalue,
            "rho": rho_vals,
            "psi_raw": psi_vals,
            "psi_norm": norm_psi,
            "r": r_vals
        }
    
    
    # ****** END OF PART 7: Solve and Cache All Wavefunctions Once ******
    
    # ****** PART 8: Solve Schrodinger Equation, Normalise, and Plot ******
    
    plt.figure()
    for i in wavefunction_data:
        rho_values = wavefunction_data[i]["rho"]
        normalized_psi = wavefunction_data[i]["psi_norm"]
        plt.plot(rho_values, normalized_psi, label=f'n={i}')
    
    plt.xlabel('ρ (dimensionless)')
    plt.ylabel('Normalized Ψ(ρ)')
    plt.title(f'Normalized Wavefunctions: {system_name}, ℓ = {ell}')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # ****** END OF PART 8: Solve Schrodinger Equation, Normalise, and Plot ******
    
    # ****** PART 9: |Ψ(0)|² and mean radius ******
    
    print("Now Calculating |Ψ(0)|² and r₀ (mean radius) for each n")
    
    for i in wavefunction_data:
        eigenvalue = wavefunction_data[i]["epsilon"]
        rho_values = wavefunction_data[i]["rho"]
        normalized_psi = wavefunction_data[i]["psi_norm"]
        r_values = wavefunction_data[i]["r"]
    
        print(f"\nState n={i} with ε={eigenvalue:.6f}")
    
        # 1. |Ψ(0)|² in ρ and r
        psi0_rho = normalized_psi[0]
        psi0_squared_rho = abs(psi0_rho)**2
        psi0_squared_r = psi0_squared_rho / (hbar_ev / (m * c))
    
        print(f"|Ψ(0)|² (in ρ): {psi0_squared_rho:.5e}")
        print(f"|Ψ(0)|² (in r): {psi0_squared_r:.5e} m⁻³")
    
        # 2. Mean radius via expectation value
        mean_r_rho = np.trapz(rho_values**3 * normalized_psi**2, rho_values)
        mean_r = mean_r_rho * (hbar_ev / (m * c))
        print(f"⟨r⟩ (expectation value): {mean_r:.5e} m")
    
        # 3. Width at half height
        half_max = max(abs(normalized_psi)) / 2
        indices = np.where(abs(normalized_psi) >= half_max)[0]
    
        if len(indices) >= 2:
            width_half_height_rho = rho_values[indices[-1]] - rho_values[indices[0]]
            width_half_height_r = width_half_height_rho * (hbar_ev / (m * c))
            print(f"r₀ (width at half max): {width_half_height_r:.5e} m")
        else:
            print("Could not determine r₀ (half-width)")
    
    # ****** END OF PART 9: |Ψ(0)|² and mean radius ******
    
    # ****** PART 10: Conversion to Standard Units ******
    
    print("\n Wavefunctions and Energies in Dimensional r ")
    
    plt.figure()
    
    for i in wavefunction_data:
        epsilon = wavefunction_data[i]["epsilon"]
        r_values = wavefunction_data[i]["r"]
        normalized_psi = wavefunction_data[i]["psi_norm"]
    
        energy_ev = epsilon * m # Again, mass is already in eV, no need to multiply by c^2
        energy_mev = energy_ev / 1e6
    
        plt.plot(r_values * 1e15, normalized_psi, label=f'n={i} ({energy_mev:.1f} MeV)')
        print(f"State n={i}: E = {energy_mev:.2f} MeV")
    
    plt.xlabel('r (fm)')
    plt.ylabel('Normalized Ψ(r)')
    plt.title(f'Wavefunctions in Dimensional r: {system_name}, ℓ = {ell}')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # ****** END OF PART 10: Conversion to Standard Units ******
    
    # Store data for combined heatmap plotting
    key = f"{system_name}, ℓ={ell}"
    all_wavefunction_data[key] = wavefunction_data


# ****** END OF CONTROL LOOP FOR SYSTEMS ******

# ****** PART 11: Heatmap Graph ******

# Plot a 2×2 grid of wavefunction heatmaps for all systems
def plot_all_wavefunction_heatmaps(data_dict):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for ax, (label, wavefunction_data) in zip(axes, data_dict.items()):
        # Interpolate wavefunctions onto a common r grid
        r_values = wavefunction_data[0]["r"]
        grid_r = np.linspace(r_values.min(), r_values.max(), 400)
        heatmap = []

        for n in range(5):
            psi = wavefunction_data[n]["psi_norm"]
            r_n = wavefunction_data[n]["r"]
            psi_interp = np.interp(grid_r, r_n, psi)
            heatmap.append(psi_interp)

        heatmap = np.array(heatmap)

        im = ax.imshow(heatmap, extent=[grid_r[0]*1e15, grid_r[-1]*1e15, 4, 0],
                       aspect='auto', cmap='viridis')
        ax.set_title(label, fontsize=10)
        ax.set_xlabel('r (fm)', fontsize=9)
        ax.set_ylabel('State index n', fontsize=9)
        ax.tick_params(labelsize=8)

    # Add shared colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label='Ψ(r)')

    plt.suptitle("Wavefunction Heatmaps for Charmonium & Bottomonium (S and P waves)", fontsize=14, weight='bold')
    plt.tight_layout(rect=[0, 0, 0.88, 0.95])
    plt.show()

# Call the function to draw the combined heatmap
plot_all_wavefunction_heatmaps(all_wavefunction_data)

# ****** END OF PART 11: Heatwave Graph ******

'''
END OF SCRIPT
'''