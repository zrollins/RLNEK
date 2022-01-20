###This is an example of RLNEK Module 0B prompts & inputs###
###This is meant to ensure RLNEK was installed correctly (DISCLAIMER: not actual receptor-ligand kinetics)###
###Please check your outputs to correspond to outputs in the outputs folder of Module_0B_test###
###Note, when NOT for testing, always copy Module_0B.py to the directory where Trackmate files are located###

###start by changing the operating directory in 'Module_0B.py' --> os.chdir('/User/RLNEK_tests/Module_0B_test')###
###Next, run script & follow prompts with described inputs below for testing###
###prompts precede ":" and user inputs follow ":"###

Enter fluid viscosity (dyne-s/cm²): 6.92e-3

Enter cell/sphere radius (μm): 7.85

Enter critical distance (μm): 0.05

Enter receptor-ligand bond length (nm): 36

Enter flow chamber height (μm): 75

Enter flow chamber width (μm): 800

Enter CCD FPS: 8

Enter minimum displacement (enter for 0.5 μm): D_min (μm) = 0.5

Enter "y" if ligand site density was characterized; otherwise, enter "n": n

Enter ligand coating concentration (μg/mL): 10

Enter flow rate (μL/hr): 5

For flow rate = 5.000000 (μL/hr), enter name of "spots" file(s) (from Trackmate) for non-specific ligand: 10ugmL_receptor_ns_ligand_5uLhr_8FPS_Trial_1_spots, 10ugmL_receptor_ns_ligand_5uLhr_8FPS_Trial_2_spots, 10ugmL_receptor_ns_ligand_5uLhr_8FPS_Trial_3_spots

For flow rate = 5.000000 (μL/hr), enter name of "spots" file(s) (from Trackmate) for specific ligand: 10ugmL_receptor_s_ligand_5uLhr_8FPS_Trial_1_spots, 10ugmL_receptor_s_ligand_5uLhr_8FPS_Trial_2_spots, 10ugmL_receptor_s_ligand_5uLhr_8FPS_Trial_3_spots

Non-specific bond lifetime (seconds): 35.714470 ± 2.470085
Specific bond lifetime (seconds): 55.095972 ± 2.159925
p=0.000000
Specific bond lifetimes are significantly GREATER than non-specific bond lifetimes. Your data looks reasonable!
