# Heavy Meson States — Kevin McCarthy, with Daniel Farrell

## Abstract
This project investigated the nature of Charm and Beauty (Bottom) Quark meson states, with particular emphasis on their respective potentials.
This was done by numerically solving the two body Schrödinger Equation eigenvalue problem with the Cornell potential, giving us a method of evaluating wavefunctions, 
and thus, mean radii and decay constants.
When comparing with data from the Particle Data Group archive, 
the larger mass of the Beauty quark was shown to allow greater computational accuracy in the non-relativistic approximation than that of Charm.
Further improvements to our theoretical model are planned for a more refined analysis in the future.


## Repository Contents
| File | Description |
|------|-------------|
| `Charm_and_Beauty.py` | Main analysis script — solves Schrodinger equations and calculates wavefunctions |
| `Data_Comparisson_Tables.py` | Uses data from Particle Data Group to compare against numerical results |
| `Energy_Difference_Comparissons.py` | Plots comparison graphs |
| `The_Potential_of_Charm_and_Beauty.pdf` | Full write-up of methodology and results |

## Tools
- Python 3.13
- Libraries used: NumPy, SciPy, Matplotlib, Pandas


## Key Results
This project presented a comprehensive numerical study of heavy quarkonium bound states using the radial Schrödinger Equation with a Cornell Potential, implemented through Python code. Key outputs included:  
•	Eigenenergies for charm and bottom quark systems,  
•	Normalized wavefunctions,  
•	Radial expectation values,  
•	Wavefunction values at the origin and spatial width,  
•	Comparisons against PDG experimental data.  
While the code correctly captured qualitative trends — such as energy level ordering, wavefunction behaviour under angular momentum, and size reduction in heavier systems — it fell short of quantitative accuracy. The discrepancies can be traced to several sources which include:  
•	Non-relativistic treatment and neglect of spin effects,  
•	Simplistic numerical integration and step size rigidity,  
•	Missing QCD radiative corrections in decay-related quantities.  
 Nonetheless, the structure of the code, including modular components for solving the Schrödinger equation, finding eigenvalues via bracket-and-bisection methods, and normalizing wavefunctions, offers a solid foundation for further development.  

## Potential Improvements
Upon reflection, it is clear to us that there are many plausible methodologies which could be used to improve our report process. 
These are mostly theoretical in nature, founded upon already known ideas surrounding our topic, physics of subatomic particles. 
Firstly, we would definitely introduce relativistic corrections to our data via Bret-Fermi or Dirac Formulism. 
Despite the solution of a non-relativistic Schrödinger Equation giving great insight into our energy spectra and other collected data, 
relativistic theory would have added to our overall accuracy greatly. 
Additionally, we would consider fixing our potential parameters to that of known experimental energy levels. 
This would give us greater sources of data to compare and contrast our script data to, and project somewhat, an ‘end goal’ for our work. 
That is, if our work garnered similar values of measurement, we would surely be quite accurate in our investigation. 
In our attempt to get higher precision, electing to utilise adaptive integration and/or spectral methods would have greatly helped us. 
Perhaps we are limited by our experience and knowledge with respect to computational physics to truly explore the possible accuracy we are capable of in this particular investigation. 
The incorporation of spin-dependent interactions for full-spectral splitting would have been of interest. 
Spin is a particularly important idea in particle physics and physically would have been crucial in raising our overall accuracy. 
Finally, calibrating the wavefunction normalisation to physical units using experimental observables was another conclusion we arrived at. 
Doing so would produce more physically relevant solutions. 
With these enhancements, the model we have produced could evolve into a powerful tool for exploring meson spectra, decay processes, 
and beyond Standard Model scenarios involving bound-state dynamics. 


## How to Run
```bash
pip install numpy scipy matplotlib
python Charm_and_Beauty.py
```

## Author
Kevin McCarthy | mccarthykevin012@gmail.com |
