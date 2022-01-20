Corresponding Manuscript:
"Receptor-Ligand Non-Equilibrium Kinetics (RLNEK) 1.0: An Integrated Trackmate Laminar Flow Chamber Analysis"

Trackmate Dependencies:
>=ImageJ 2.3.0/1.53f51
>=Java 1.8.0_301 [64-bit]
>=Trackmate v7.4.0

RLNEK Dependencies:
>= Python 3.7
>= pandas 1.3.5
>= numpy 1.22.0
>= scipy 1.7.3
>= matplotlib 3.5.1

RLNEK Modules:
Module 0A - convert flow rate to applied tether/tensile force
Module 0B - determine non-specific bond lifetime
Module 0C - determine if coating concentration/site density range meets single-molecule criteria
Module 1 - determine Nb/NT & koff @ various flow rates, coating concentrations
NT Optimization - optimize NT values for given flow rate across various site densities (precursor to Module 2)
Module 2 - determine Nb/NT, koff, k+, & kin @ various flow rates, site densities


Tips:
1. Please read manuscript and corresponding User Guide to understand experimental parameters and data necessary to perform calculations.
2. In previous versions of Trackmate, the files are outputted in different formats. Ensure that Trackmate "tracks" and "spots" files are formatted correctly (e.g., in the "spots" file, particleIDS are listed in numerical order of trackIDs & numerical order of frames). Use provided test files as a reference.
3. All modules can handle replicates for a given condition. Please see test folders for implementation.
4. Collect Receptor-Ligand Non-Equilibrium Kinetics.

Appendix:
Nb = number of cells/spheres bound
NT = total number of cells/spheres to cross field of view
Nb/NT = capture efficiency
koff = force-dependent dissociation rate (1/ bond lifetime)
k+ = effective on rate of cell/spheres (includes diffusive & kinetic timescales)
kin = intrinsic on rate of receptor-ligand pair (assumed constant)