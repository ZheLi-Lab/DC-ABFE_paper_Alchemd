# Alchemd: 
A versatile and easy to use software for FEP calculation.

## Installation
1. Install Dependencies:
- Python 3.9 or higher
- cudatoolkit 11.0
- openmm 7.7
- openmmtools 0.21.5
- plumed 2.9.2
- openmm-plumed 1.0
- tqdm 4.67.1
- parmed 
- alive-progress 2.4.1
- matplotlib 3.6.0
```sh
conda create --name alchemd_fep python=3.9
conda activate alchemd_fep
conda install openmm=7.7 openmmtools=0.21.5 cudatoolkit=11.0 -c conda-forge
conda install plumed openmm-plumed tqdm -c conda-forge
conda install parmed -c conda-forge
pip install alive-progress==2.4.1
pip install matplotlib==3.6.0
```

2. Install Alchemd:
```sh
git clone https://github.com/***
cd Alchemd
pip install -e .
```

## Basic Usage and Some Basic Control files
1. Basic Usage:
```sh
run_alchemd -p protein.prmtop -c protein.rst7 -i input.txt > mdout 2>&1
```
where protein.prmtop is the amber style topology file, protein.rst7 is the amber style coordinate file (could be protein.prmcrd files, but the suffix need to be changed to .rst7), and input.txt is the input control file.

2. Some Basic Control files:
- lambdas.json:
```json
{
    "res_lambda": 
    {
        "lambda_restraints": [0.0, 1.0], 
        "lambda_electrostatics": [1.0, 1.0], 
        "lambda_sterics": [1.0, 1.0]
    }, 
    "cal_mbar_lambda": 
    {
        "lambda_restraints": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 
        "lambda_electrostatics": [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0], 
        "lambda_sterics": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    }
}
```
The lambdas.json file contains the lambda values for the defined states (whose names are defined as the "lambdas_json" option in the input.txt file). There may be multiple lambda dictionaries in the file. The user can use the "lambdas_group" option in the input.txt file to specify which lambda dictionary to use.
**Note:** For OpenMM, when the lambda_electrostatics and lambda_sterics values are set to 1.0, the ligand will be fully interacting with the environment with electrostatics and vdW interactions. When the lambda_electrostatics and lambda_sterics values are set to 0.0, the interaction between the ligand and environment will be fully turned off. However, for molecular dynamics simulation software like Gromacs, the definition of interaction states corresponding to lambda_electrostatics and lambda_sterics values of 1 and 0 is exactly opposite to that of OpenMM. To maintain consistency with Gromacs and similar software, the output files (state_s*.csv) generated by Alchemd adopt the same definition as Gromacs regarding the relationship between λ values and interaction states. However, in the lambdas.json file that serves as input for OpenMM, the OpenMM-style definition of interaction states corresponding to λ=1 or 0 is retained.

- input.txt:
-- Scenario 1: Run the preliminary MD simulation and use the generated sampling data and the RED-E based restraint selection method to define the Boresch-like restraints.(In this scenario, the ligand residue name should be renamed to "MOL" for default option)
```
[normalmd]
normalmd=False # do not run normal md

[alchemical]
alchemical = False # do not run alchemical md

[restraint]
restraint = True
temperature = 300.0 
timestep = 2 # in fs
ligand_resname = MOL 
iflog_restraint_detail = True # if True, save the detailed restraint information in the log file (all_lst)
# increase temperature 50 times within 25000 timesteps (in this case, each time increase 6 Kelvin)
heat_nsteps = 25000 # which is 50 ps
heat_iterations = 50 # use 50 iterations to increase the temperature to 300 K
density_nsteps = 25000 # which is 50 ps
npt_nsteps = 1000000 # production with plumed, which is 2 ns
f_npt_state = ./npt_final_state.xml
save_traj = False
reportstate_freq = 1000
savetraj_freq = 50000
f_plumed_input = plumed.dat # plumed input file, will be generated automatically
f_plumed_output = Colvar # plumed output file, will be generated during the simulation
lambdas_json = lambdas.json # the lambdas.json file that contains the lambda values for the lambda_restraints states
lambdas_group = res_lambda # the group name in the lambdas.json file that contains the lambda values for the lambda_restraints states
fake_state_xml = fake_state.xml # for skipping the running of the state_s0, if not state_s0.xml, it will not work.
first_state_csv = npt_state.csv # according to the lambdas.json file and the respective lambdas_group, to generate the potential energy for each lambda_restraints state
plumed_freq = 100 # plumed record frequency
crd = None # assigned by the openmm-FEP-run.py -c option, but you can give the coordinate file of the structure that you like for the generation of candidate groups of restraint
top = None # assigned by the openmm-FEP-run.py -p option, but you can give the topology file of the structure that you like for the generation of candidate groups of restraint
res_sele_strategy = HB_pair|HB_mainchain # could be 'lig_shape|HB_pair|HB_mainchain|Huggins'
opt_cost_name = RED_E_cost # could be 'RED_E_cost' or 'dG_forward' or False
if_mean = False  # if_mean is True and the opt_cost_name is False, use the mean values of six parameters as the eq values for the harmonic restraint
if_init_pose = False # if_mean is True and the opt_cost_name is False , use the initial values of six parameters from the input pose as the eq values for the harmonic restraint
preliminary_md_inital_state = None
preliminary_min_and_heat = True # if True, will run the minimization and heating process before the production MD simulation
f_restraint = ./res_databystd.csv # the file logging the six-degree-of-freedom restraint pararmeters
```

-- Scenario 2: Run the alchemical MD simulation for classical FEP-ABFE calculation with the defined Boresch-like restraints.
```
[normalmd]
normalmd=False # do not run normal md

[alchemical]
alchemical = True # run alchemical md
lambda_json = lambdas.json # the lambdas.json file that contains the lambda values for the alchemical transformation
lambda_group = cal_mbar_lambda # the group name in the lambdas.json file that contains the lambda values for the alchemical transformation
temperature = 300.0 # in K
timestep = 2 # in fs
nsteps = 100 # the number of steps for one iteration of the alchemical transformation, which is the sampling frequency of the potential energy results
niterations = 10000 # the number of iterations for the alchemical transformation, in this case, 10000 iterations will give 10000 sampling points, the simulation time will be 10000 * 100 * 2 fs = 200000 fs
save_traj = False
input_state = ../npt_final_state.xml # the initial state of the alchemical transformation
cal_adj_ene_num = 5 # To set the number of the adjacent windows of current windows to calculate the internal energy, eg. cal_adj_ene_num=5, normally will be 11 windows energy calculated, if cal_adj_ene_num=all, calculate all windows energy.
kbond = 10 # kbond is the spring constant for the harmonic bonds, in kcal/mol/A^2
ktheta1 = 10 # ktheta1 is the spring constant for the harmonic angles, in kcal/mol/rad^2
ktheta2 = 10 # ktheta2 is the spring constant for the harmonic dihedrals, in kcal/mol/rad^2
kphi1 = 10 # kphi1 is the spring constant for the harmonic dihedrals, in kcal/mol/rad^2
kphi2 = 10 # kphi2 is the spring constant for the harmonic dihedrals, in kcal/mol/rad^2
kphi3 = 10 # kpsi is the spring constant for the harmonic dihedrals, in kcal/mol/rad^2

[restraint]
restraint = True 
f_npt_state = ../npt_final_state.xml # the final state of the preliminary MD simulation
f_restraint = ../res_databystd.csv # the file logging the six-degree-of-freedom restraint pararmeters
```

-- Scenario 3: Run the alchemical MD simulation for Divide-and-Conquer-ABFE(DC-ABFE) calculation with the defined Boresch-like restraints.
```
[normalmd]
normalmd=False # do not run normal md

[alchemical]
alchemical = True # run alchemical md
lambda_json = lambdas.json # the lambdas.json file that contains the lambda values for the alchemical transformation
lambda_group = cal_mbar_lambda # the group name in the lambdas.json file that contains the lambda values for the alchemical transformation
temperature = 300.0 # in K
timestep = 2 # in fs
nsteps = 100 # the number of steps for one iteration of the alchemical transformation, which is the sampling frequency of the potential energy results
niterations = 10000 # the number of iterations for the alchemical transformation, in this case, 10000 iterations will give 10000 sampling points, the simulation time will be 10000 * 100 * 2 fs = 200000 fs
save_traj = False
input_state = ../npt_final_state.xml # the initial state of the alchemical transformation
cal_adj_ene_num = 5 # To set the number of the adjacent windows of current windows to calculate the internal energy, eg. cal_adj_ene_num=5, normally will be 11 windows energy calculated, if cal_adj_ene_num=all, calculate all windows energy.
kbond = 10 # kbond is the spring constant for the harmonic bonds, in kcal/mol/A^2
ktheta1 = 10 # ktheta1 is the spring constant for the harmonic angles, in kcal/mol/rad^2
ktheta2 = 10 # ktheta2 is the spring constant for the harmonic dihedrals, in kcal/mol/rad^2
kphi1 = 10 # kphi1 is the spring constant for the harmonic dihedrals, in kcal/mol/rad^2
kphi2 = 10 # kphi2 is the spring constant for the harmonic dihedrals, in kcal/mol/rad^2
kphi3 = 10 # kpsi is the spring constant for the harmonic dihedrals, in kcal/mol/rad^2
pdbx = mol.pdbx # the pdbx file of the ligand
current_group_nb  = 1 # the index (starting from 1) of the current group of atoms, whose vdW interaction with the environment to be decoupled.
current_group_chg = 0 # beta-test feature, the index (starting from 1) of the columns of the charge information, which will be set to the initial charge of the ligand. The next column is the target charge of the ligand. For the decoupling the only vdW interaction by groups, the current_group_chg could be set to 0, with the lambda_electrostatics being 0.0. 

[restraint]
restraint = True 
f_npt_state = ../npt_final_state.xml # the final state of the preliminary MD simulation
f_restraint = ../res_databystd.csv # the file logging the six-degree-of-freedom restraint pararmeters
```

- pdbx file:
A example pdbx file could be:
```
HETATM 2122  C1  MOL   128      -7.553   9.633  -4.325   0.000000   0.000000 12
HETATM 2123  C2  MOL   128      -8.191  10.360  -3.288   0.000000   0.000000 11
HETATM 2124  H1  MOL   128      -8.512  11.377  -3.459   0.000000   0.000000 11
HETATM 2125  C3  MOL   128      -8.417   9.785  -2.022   0.000000   0.000000 13
HETATM 2126 Cl1  MOL   128      -9.198  10.678  -0.773   0.000000   0.000000 1
HETATM 2127  C4  MOL   128      -8.017   8.465  -1.772   0.000000   0.000000 14
HETATM 2128  H2  MOL   128      -8.187   8.018  -0.802   0.000000   0.000000 14
HETATM 2129  C5  MOL   128      -7.426   7.710  -2.794   0.000000   0.000000 17
HETATM 2130  H3  MOL   128      -7.174   6.682  -2.586   0.000000   0.000000 17
HETATM 2131  C6  MOL   128      -7.178   8.278  -4.064   0.000000   0.000000 18
HETATM 2132  N1  MOL   128      -6.639   7.444  -5.080   0.000000   0.000000 19
HETATM 2133  C7  MOL   128      -5.668   6.512  -5.022   0.000000   0.000000 20
HETATM 2134  C8  MOL   128      -4.834   6.207  -3.817   0.000000   0.000000 16
HETATM 2135  H4  MOL   128      -3.851   5.847  -4.124   0.000000   0.000000 16
HETATM 2136  H5  MOL   128      -4.690   7.096  -3.204   0.000000   0.000000 16
HETATM 2137  H6  MOL   128      -5.309   5.431  -3.217   0.000000   0.000000 16
HETATM 2138  N2  MOL   128      -5.517   5.897  -6.202   0.000000   0.000000 20
HETATM 2139  N3  MOL   128      -6.419   6.495  -7.075   0.000000   0.000000 20
HETATM 2140  C9  MOL   128      -7.054   7.422  -6.352   0.000000   0.000000 15
HETATM 2141 C10  MOL   128      -8.081   8.387  -6.851   0.000000   0.000000 7
HETATM 2142  H7  MOL   128      -8.313   8.150  -7.889   0.000000   0.000000 7
HETATM 2143  H8  MOL   128      -9.008   8.268  -6.289   0.000000   0.000000 7
HETATM 2144  N4  MOL   128      -7.575   9.745  -6.765   0.000000   0.000000 6
HETATM 2145 C11  MOL   128      -7.323  10.318  -5.647   0.000000   0.000000 10
HETATM 2146 C12  MOL   128      -6.750  11.706  -5.711   0.000000   0.000000 9
HETATM 2147 C13  MOL   128      -6.867  12.475  -6.897   0.000000   0.000000 3
HETATM 2148  H9  MOL   128      -7.411  12.097  -7.751   0.000000   0.000000 3
HETATM 2149 C14  MOL   128      -6.256  13.740  -7.001   0.000000   0.000000 2
HETATM 2150 H10  MOL   128      -6.343  14.308  -7.916   0.000000   0.000000 2
HETATM 2151 C15  MOL   128      -5.518  14.259  -5.922   0.000000   0.000000 4
HETATM 2152 H11  MOL   128      -5.039  15.223  -6.008   0.000000   0.000000 4
HETATM 2153 C16  MOL   128      -5.402  13.515  -4.735   0.000000   0.000000 5
HETATM 2154 H12  MOL   128      -4.827  13.906  -3.908   0.000000   0.000000 5
HETATM 2155 C17  MOL   128      -6.019  12.254  -4.627   0.000000   0.000000 8
HETATM 2156 H13  MOL   128      -5.883  11.699  -3.712   0.000000   0.000000 8
```
**Note:** 
Although pdbx files are similar to pdb files, pdbx files are parsed using spaces as delimiters. Each atom occupies one line, and each line primarily contains the following information: atom number (which must correspond to the atom numbers in the prmtop and rst7 files being used), atom name, residue name, residue number, and x, y, z coordinates. (This information largely corresponds to what's available in pdb files).
However, pdbx files contain two additional types of information. One is the vdW annihilation group that each atom belongs to (corresponding to current_group_nb), which is specified in the last column of each space-delimited line. The other additional information consists of all floating-point numbers between the z-coordinate and the last element, which are parsed and combined into a list. We refer to this list as the charge information for each atom. We can specify the charge information for the initial state using current_group_chg (starting from 1), and use the charge information from column (current_group_chg+1) as the target state's charge information.
The control functionality for charge information is currently in beta development.