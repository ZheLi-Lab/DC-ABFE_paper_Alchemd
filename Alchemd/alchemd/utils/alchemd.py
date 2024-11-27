try:
    import openmm
    from openmm import unit
    from openmm import app
except ImportError:  # Openmm < 7.6
    from simtk import openmm
    from simtk import unit
    from simtk.openmm import app

import copy,sys
import sys,os
import numpy as np
import mdtraj
from openmmtools import integrators, states, cache
from openmmtools.forcefactories import restrain_atoms_by_dsl
#from openmmtools import alchemy
from . import alchemy_by_group as alchemy
from .Restraints_Select.restraints import Boresch, RestraintParameterError, RestraintState
from .Restraints_Select.restraints2 import Boresch2, RestraintState2 
from .tools.tools import Timer
import pandas as pd
from .pdbx.pdbx_parser import PdbxParser, InteractionGroups


def get_adjacent_numbers(n, threshold, adjacent_num=5, num_pre_adjacent=1) -> list:
    """
    Generate a list of numbers adjacent to a given value within specified bounds, 
    with support for multi-number state blocks.

    Parameters
    ----------
    n : int
        The central number around which to generate adjacent values
    threshold : int
        The upper bound limit for generated numbers
    adjacent_num : int, optional
        Number of adjacent states to generate on each side (default=5)
    num_pre_adjacent : int, optional
        Number of consecutive numbers that constitute one state block (default=1).
        For example:
        - num_pre_adjacent=1: Each number represents one state
        - num_pre_adjacent=2: Every two consecutive numbers represent one state

    Returns
    -------
    list
        List of integers containing n and adjacent numbers within bounds [0, threshold],
        organized according to the num_pre_adjacent parameter

    Examples
    --------
    >>> # With num_pre_adjacent=1 (default)
    >>> get_adjacent_numbers(5, 10, 2)
    [3, 4, 5, 6, 7]
    
    >>> # With num_pre_adjacent=2 (two numbers per state)
    >>> get_adjacent_numbers(5, 10, 2, 2)
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    
    Notes
    -----
    When num_pre_adjacent > 1, the function generates numbers for complete state blocks,
    which may result in more numbers than with num_pre_adjacent=1
    """
    result = []
    for i in range(max(0, n - adjacent_num * num_pre_adjacent),
                   min(threshold, (n + num_pre_adjacent - 1) + adjacent_num * num_pre_adjacent) + 1):
        result.append(i)
    return result

def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.3+
    """
    Create an iterator with a progress bar for tracking iteration progress.

    Parameters
    ----------
    it : iterable
        The iterable object to track progress for
    prefix : str, optional
        Text to display before the progress bar (default="")
    size : int, optional
        Width of the progress bar in characters (default=60)
    out : file object, optional
        Output stream for the progress bar (default=sys.stdout)

    Yields
    ------
    item
        Items from the input iterable

    Notes
    -----
    Displays a progress bar with the format:
    prefix [####....] current/total

    Examples
    --------
    >>> for i in progressbar(range(100)):
    ...     time.sleep(0.1)
    [##########] 100/100
    """
    count = len(it)
    def show(j):
        x = int(size*j/count)
        print("{}[{}{}] {}/{}".format(prefix, "#"*x, "."*(size-x), j, count), 
                end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)

class RestraintParam():
    """
    Store and manage parameters for Boresch-style geometric restraints between receptor and ligand atoms.

    This class handles the parameters needed for applying positional restraints in
    alchemical free energy calculations, including distances, angles, dihedrals and
    their associated force constants.

    Parameters
    ----------
    rec_atoms : list
        List of 3 receptor atom indices [rec1, rec2, rec3] for restraint definition, the atom number starts from 0.
    lig_atoms : list
        List of 3 ligand atom indices [lig1, lig2, lig3] for restraint definition, the atom number starts from 0.
    r : float
        Distance between rec3 and lig1 in Angstroms
    theta1 : float
        Angle between rec2-rec3-lig1 in radians
    theta2 : float
        Angle between rec3-lig1-lig2 in radians
    phi1 : float
        Dihedral angle between rec1-rec2-rec3-lig1 in radians
    phi2 : float
        Dihedral angle between rec2-rec3-lig1-lig2 in radians
    phi3 : float
        Dihedral angle between rec3-lig1-lig2-lig3 in radians
    kbond : float, optional
        Force constant for distance restraint in kcal/mol/Å² (default=10)
    ktheta1 : float, optional
        Force constant for theta1 angle in kcal/mol/rad² (default=10)
    ktheta2 : float, optional
        Force constant for theta2 angle in kcal/mol/rad² (default=10)
    kphi1 : float, optional
        Force constant for phi1 dihedral in kcal/mol/rad² (default=10)
    kphi2 : float, optional
        Force constant for phi2 dihedral in kcal/mol/rad² (default=10)
    kphi3 : float, optional
        Force constant for phi3 dihedral in kcal/mol/rad² (default=10)

    Notes
    -----
    For single atom restraints, ktheta2, kphi2, and kphi3 should be set to zero.
    All angles should be provided in radians.

    Reference
    --------
    1. S. Boresch, F. Tettinger, M. Leitgeb and M. Karplus, Absolute Binding Free Energies: A Quantitative Approach for Their Calculation, J. Phys. Chem. B, 2003, 107(35), 9535–9551, DOI: 10.1021/jp0217839.

    Examples
    --------
    >>> # Define restraints between receptor atoms [0,1,2] and ligand atoms [3,4,5]
    >>> restraints = RestraintParam([0,1,2], [3,4,5], r=5.0, theta1=1.57,
    ...                            theta2=1.57, phi1=0.0, phi2=0.0, phi3=0.0)
    """
    def __init__(self, rec_atoms, lig_atoms, r, theta1, theta2, phi1, phi2, phi3, kbond=10, ktheta1=10, ktheta2=10, kphi1=10, kphi2=10, kphi3=10):
        self.rec_atoms = rec_atoms # [rec1, rec2, rec3]
        self.lig_atoms = lig_atoms # [lig1, lig2, lig3]
        self.r = r # distance between rec3 and lig1, unit in angstrom
        self.theta1 = theta1 # angle between rec2, rec3, lig1, unit in radians 
        self.theta2 = theta2 # angle between rec3, lig1, lig2  
        self.phi1 = phi1 # dihedral between rec1, rec2, rec3, lig1
        self.phi2 = phi2 # dihedral between rec2, rec3, lig1, lig2 
        self.phi3 = phi3 # dihedral between rec3, lig1, lig2, lig3 
        self.kbond = kbond # unit in kcal/mol/A^2
        self.ktheta1 = ktheta1 # unit in kcal/mol/rad^2
        self.ktheta2 = ktheta2 # unit in kcal/mol/rad^2 for single atom restraint, this term to be zero
        self.kphi1 = kphi1 # unit in kcal/mol/rad^2
        self.kphi2 = kphi2 # unit in kcal/mol/rad^2 for single atom restraint, this term to be zero
        self.kphi3 = kphi3 # unit in kcal/mol/rad^2 for single atom restraint, this term to be zero
    def __repr__(self):
        return f'rec_atoms:{self.rec_atoms}, lig_atoms:{self.lig_atoms}, r:{self.r} A, theta1:{self.theta1} rad, theta2:{self.theta2} rad, phi1:{self.phi1} rad, phi2:{self.phi2} rad, phi3:{self.phi3} rad, kbond:{self.kbond} kcal/mol/A^2, K_theta1:{self.ktheta1} kcal/mol/rad^2, self.K_theta2:{self.ktheta2} kcal/mol/rad^2, self.K_phi1:{self.kphi1} kcal/mol/rad^2, self.K_phi2:{self.kphi2} kcal/mol/rad^2, self.K_phi3:{self.kphi3} kcal/mol/rad^2'

class MDClass(app.Simulation):
    """
    A molecular dynamics simulation class extending OpenMM's app.Simulation.

    This class provides additional functionality for molecular dynamics simulations,
    including state handling, PLUMED integration.

    Parameters
    ----------
    system_container : object
        Container object that must have attributes:
        - topology : OpenMM Topology object
        - system : OpenMM System object

    Attributes
    ----------
    topology : openmm.app.Topology
        Molecular topology of the system
    system : openmm.System
        OpenMM system containing forces and constraints
    reporters : list
        List of reporter objects for trajectory/data output
    _usesPBC : bool
        Flag indicating if periodic boundary conditions are used

    Methods
    -------
    set_state_xpv(state)
        Set positions, velocities and box vectors from a State object
    useplumed(plumed_filename)
        Add PLUMED forces defined in input file to the system

    Properties
    ----------
    final_state : openmm.State
        Current state with positions, velocities and box vectors
    sampler_state : openmmtools.states.SamplerState
        Current state formatted for enhanced sampling frameworks

    Examples
    --------
    >>> md = MDClass(system_container)
    >>> md.useplumed("plumed.dat")
    >>> md.step(1000)
    """
    def __init__(self, system_container):
        self.topology = system_container.topology 
        self.system = system_container.system 
        self.reporters = []
        self._usesPBC = True

    def set_state_xpv(self,state):
        """
        Set context state from an OpenMM State object.

        Parameters
        ----------
        state : openmm.State
            State object containing positions, velocities and box vectors
        """
        #set context positions, velocities, box vectors based on state
        positions = state.getPositions()
        velocities = state.getVelocities()
        boxvectors = state.getPeriodicBoxVectors()
        self.context.setPositions(positions)
        self.context.setVelocities(velocities)
        self.context.setPeriodicBoxVectors(*boxvectors)

    def useplumed(self,plumed_filename):
        """
        Add PLUMED forces to the system from an input file.

        Parameters
        ----------
        plumed_filename : str
            Path to PLUMED input file defining biasing forces/CVs

        Notes
        -----
        This method modifies the system by adding a PlumedForce
        """
        print("Using plumed file %s ..." % (plumed_filename))
        try:
            from openmmplumed import PlumedForce
        except:
            print("Error: Failed to import openmmplumed. Please make sure it is installed correctly.")
            sys.exit(1)
        plumed_infile = open(plumed_filename,'r')
        try:
            plumed_script = plumed_infile.read()
        finally:
            plumed_infile.close()
        self.plumedforce = PlumedForce(plumed_script)
        self.system_container.system.addForce(self.plumedforce)

    @property
    def final_state(self):
        return self.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)

    @property
    def sampler_state(self):
        return states.SamplerState.from_context(self.context)


class NormalMD(MDClass):
    """
    Standard molecular dynamics simulation with heating and equilibration protocols.

    A class for running molecular dynamics simulations with support for minimization,
    heating, and density equilibration, inheriting from MDClass.

    Parameters
    ----------
    system_container : object
        Container with system topology, positions and forces
    timestep : openmm.unit.Quantity, optional
        Integration timestep (default=4*femtoseconds)
    temperature : openmm.unit.Quantity, optional
        Target temperature (default=298*kelvin)
    constant_pressure : openmm.unit.Quantity, optional
        Pressure for NPT simulation (default=None for NVT)
    plumed : str, optional
        Path to PLUMED input file (default=None)
    read_state : str or openmm.State, optional
        Initial state to load, either file path or State object
    platform : openmm.Platform, optional
        OpenMM platform to use for computation
    platformProperties : dict, optional
        Platform-specific properties

    Methods
    -------
    minimize()
        Perform local energy minimization
    minimize_cascade(mode='com')
        Multi-step minimization with decreasing restraints
    minimize_cascade_lig()
        Minimization protocol specific for ligand systems
    minimize_cascade_com()
        Minimization protocol for protein-ligand complexes
    heat_cascade(nsteps=50000, heat_iterations=60)
        Gradual heating with restraints
    density_run(nsteps=10000, mode='com')
        Density equilibration under NPT
    run(nsteps, report_speed=True)
        Production MD with optional performance reporting

    Notes
    -----
    The minimize_cascade methods use a series of positional restraints that are
    gradually removed. The heat_cascade method gradually increases temperature
    while maintaining restraints.

    Examples
    --------
    >>> md = NormalMD(system_container, temperature=300*unit.kelvin)
    >>> md.minimize_cascade()
    >>> md.heat_cascade()
    >>> md.density_run(10000)
    >>> md.run(100000)
    """
    def __init__(self, system_container,
                timestep=4*unit.femtoseconds, temperature=298*unit.kelvin, 
                constant_pressure=None, plumed=None, read_state=None,
                platform=None, platformProperties=None):
        self.system_container=copy.deepcopy(system_container)
        self.original_system_container=system_container
        self.read_state=read_state
        self.timestep=timestep
        self.temperature=temperature
        self.constant_pressure=constant_pressure
        super().__init__(self.system_container)
        self.topology=self.original_system_container.topology #could not use self.system_container. Could be a bug
        if plumed:
            self.plumedforce = None
            self.useplumed(plumed)
        if constant_pressure:
            self.thermodynamic_state = states.ThermodynamicState(self.system_container.system, temperature, constant_pressure)
        else:
            self.thermodynamic_state = states.ThermodynamicState(self.system_container.system, temperature)
        # self.integrator= integrators.LangevinIntegrator(timestep=self.timestep,splitting="V R R R O R R R V")
        self.integrator = openmm.openmm.LangevinMiddleIntegrator(self.temperature, 1/unit.picoseconds, self.timestep)#TO VERTIFY
        self.context=self.thermodynamic_state.create_context(self.integrator,platform=platform,platform_properties=platformProperties)
        # adding positional restraints will change thermodynamic_state. Save a copy
        self.thermodynamic_state_original=copy.deepcopy(self.thermodynamic_state)

        if not read_state:
            self.context.setPositions(self.system_container.positions)
        else:
            if isinstance(read_state, str):
                self.loadState(read_state) # change this to sampler_state
            else:
                self.set_state_xpv(read_state)

    def minimize(self):
        """Perform local energy minimization of the system."""
        openmm.LocalEnergyMinimizer.minimize(self.context)

    def minimize_cascade(self, mode='com'):
        """
        Perform staged minimization with decreasing restraints.

        Parameters
        ----------
        mode : {'com', 'lig'}, optional
            Minimization mode - 'com' for protein-ligand complex (default),
            'lig' for ligand-only system
        """
        #minimize cascade,
        # mode: 'com', perform minimization for a protein-ligand complex system, default
        # mode: 'lig', perform minimization for a ligand system
        if mode == 'com':
            self.minimize_cascade_com()
        elif mode == 'lig':
            self.minimize_cascade_lig()

    def minimize_cascade_lig(self):
        """
        Perform a two-step minimization cascade specifically for ligand systems.

        The cascade consists of:
        1. Minimization with restraints on all non-water heavy atoms
        2. Minimization with no restraints

        Each step uses positional restraints (sigma=0.1Å) that are modified between
        steps. The thermodynamic state is reset between steps.

        Notes
        -----
        Restraints are implemented using the format:
        - Step 1: '(all and not water) and not (name =~ "H.*")'
        - Step 2: No restraints
        """
        print('Starting 4 steps minimization cascade')
        print('Performing 1st minimization step...')
        #adding positional restraints will change thermodynamic_state
        restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology, '(all and not water) and not (name =~ "H.*")',sigma=0.1*unit.angstrom)
        #self.thermodynamic_state.K=100.0*unit.kilocalories_per_mole/unit.angstrom**2
        self.context.reinitialize(preserveState=True)
        self.minimize()

        print('Performing 2rd minimization step...')
        #remove positional restraint forces.
        self.thermodynamic_state=copy.deepcopy(self.thermodynamic_state_original)
        self.context.reinitialize(preserveState=True)
        self.minimize()

    def minimize_cascade_com(self):
        """
        Perform a four-step minimization cascade for protein-ligand complexes.

        The cascade consists of:
        1. Minimization with restraints on all non-water heavy atoms
        2. Minimization with restraints on all non-water, non-ligand heavy atoms
        3. Minimization with restraints only on backbone heavy atoms
        4. Minimization with no restraints

        Each step uses positional restraints (sigma=0.1Å) that are modified between
        steps. The thermodynamic state is reset between steps.

        Notes
        -----
        Restraints are implemented using the format:
        - Step 1: '(all and not water) and not (name =~ "H.*")'
        - Step 2: '(all and not water and not resname MOL) and not (name =~ "H.*")'
        - Step 3: 'backbone and not (name =~ "H.*")'
        - Step 4: No restraints
        """
        print('Starting 4 steps minimization cascade')
        print('Performing 1st minimization step...')
        #adding positional restraints will change thermodynamic_state
        restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology, '(all and not water) and not (name =~ "H.*")',sigma=0.1*unit.angstrom)
        #self.thermodynamic_state.K=100.0*unit.kilocalories_per_mole/unit.angstrom**2
        self.context.reinitialize(preserveState=True)
        self.minimize()

        print('Performing 2nd minimization step...')
        #add a new type of positional restraints, first we need to restore the original thermodynamic_state
        self.thermodynamic_state=copy.deepcopy(self.thermodynamic_state_original)
        restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology, '(all and not water and not resname MOL) and not (name =~ "H.*")',sigma=0.1*unit.angstrom)
        #self.thermodynamic_state.K=100.0*unit.kilocalories_per_mole/unit.angstrom**2
        self.context.reinitialize(preserveState=True)
        self.minimize()

        print('Performing 3rd minimization step...')
        #add a new type of positional restraints, first we need to restore the original thermodynamic_state
        self.thermodynamic_state=copy.deepcopy(self.thermodynamic_state_original)
        restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology, 'backbone and not (name =~ "H.*")',sigma=0.1*unit.angstrom)
        #self.thermodynamic_state.K=100.0*unit.kilocalories_per_mole/unit.angstrom**2
        self.context.reinitialize(preserveState=True)
        self.minimize()

        print('Performing 4th minimization step...')
        #remove positional restraint forces.
        self.thermodynamic_state=copy.deepcopy(self.thermodynamic_state_original)
        self.context.reinitialize(preserveState=True)
        self.minimize()

    def heat_cascade(self,nsteps=50000,heat_iterations=60):
        """
        Gradually heat the system while maintaining restraints.

        Parameters
        ----------
        nsteps : int, optional
            Total number of MD steps (default=50000)
        heat_iterations : int, optional
            Number of temperature increments (default=60)
        """
        print('Starting heat cascade, total steps %d, heat iterations %d' % (nsteps, heat_iterations))
        self.thermodynamic_state=copy.deepcopy(self.thermodynamic_state_original)
        restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology, 'all and not water')
        self.thermodynamic_state.K=10.0*unit.kilocalories_per_mole/unit.angstrom**2
        self.context.reinitialize(preserveState=True)
        T_increment=self.temperature/heat_iterations
        for i in range(heat_iterations):
            T=(i+1)*T_increment
            self.integrator.setTemperature(T*unit.kelvin)
            self.run(int(nsteps/heat_iterations), report_speed=False)

    def density_run(self, nsteps=10000, mode='com'):
        """
        Equilibrate system density under NPT conditions.

        Parameters
        ----------
        nsteps : int, optional
            Number of MD steps (default=10000)
        mode : {'com', 'lig'}, optional
            System type - 'com' for complex, 'lig' for ligand (default='com')

        Raises
        ------
        ValueError
            If constant_pressure is not set
        """
        print('Starting density equilibration, total steps %d' % (nsteps))
        #add positional restraints
        if not self.constant_pressure:
            raise ValueError('Density run must be performed under constant pressure')
        self.thermodynamic_state=copy.deepcopy(self.thermodynamic_state_original)

        # Bug will occur when 'protein' represents more than molecules.
        # TODO: change to positional restraint may solve this problem.
        if mode == 'com':
            try:
                restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology, 'protein')
            except:
                restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology, 'resname MOL')
        elif mode == 'lig':
            restrain_atoms_by_dsl(self.thermodynamic_state, self.sampler_state, self.topology, 'resname MOL')

        self.thermodynamic_state.K=10.0*unit.kilocalories_per_mole/unit.angstrom**2
        self.context.reinitialize(preserveState=True)
        self.run(nsteps)
        #remove positional restraints
        self.thermodynamic_state.K=0.0*unit.kilocalories_per_mole/unit.angstrom**2
        self.thermodynamic_state.apply_to_context(self.context)

    def run(self,nsteps,report_speed=True):
        """
        Run production molecular dynamics.

        Parameters
        ----------
        nsteps : int
            Number of MD steps to run
        report_speed : bool, optional
            Whether to report simulation speed in ns/day (default=True)
        """
        if report_speed:
            t=Timer()
            t.start()

        self._simulate(endStep=self.currentStep+nsteps)
        
        if report_speed:
            clock_time=t.stop()
            ns_per_day=self.timestep.value_in_unit(unit=unit.nanosecond)*nsteps*86400/clock_time
            print("Simulation Speed: %4.3f ns/day" % (ns_per_day))


class AlchemMD(MDClass):
    """
    Alchemical molecular dynamics simulation class for free energy calculations.

    A class for running alchemical transformations with support for classical FEP-ABFE
    or DC-ABFE.

    Parameters
    ----------
    system_container : object
        Container with system topology, positions and forces 
    restraint_parm : RestraintParam, optional
        Parameters for restraint potentials(which is controlled by the lambda_restraints)
    timestep : openmm.unit.Quantity, optional
        Integration timestep (default=4*femtoseconds)
    temperature : openmm.unit.Quantity, optional
        Target temperature (default=298*kelvin)
    pdbx : str, optional
        Path to PDBX file specifying groups for VDW decoupling and charge 
        annihilation
    current_group_nb : int, optional
        Current non-bonded group index (default=0)
    current_group_chg : int, optional
        Current charge group index (default=0)
    constant_pressure : openmm.unit.Quantity, optional
        Pressure for NPT simulation (default=None for NVT)
    plumed : str, optional
        Path to PLUMED input file
    read_state : str or openmm.State, optional
        Initial state to load
    annihilate_electrostatics : bool, optional
        Whether to annihilate electrostatics (default=True)
    alchem_group : str, optional
        DSL selection for alchemical region (default='resname MOL')
    platform : openmm.Platform, optional
        OpenMM platform to use
    platformProperties : dict, optional
        Platform-specific properties
    another_restraint_parm : RestraintParam, optional
        Additional restraint parameters(which is controlled by the lambda_restraints2)
    **kwargs
        Additional arguments including alchem_settings

    Attributes
    ----------
    alchemical_state : AlchemicalState
        Current alchemical state
    restraint_state : RestraintState  
        Current restraint state
    u_kln_array : array
        Energy matrix for MBAR analysis
    lambdas : dict, 
        keys: list of str, like ["lambda_restraints", "lambda_restraints2", "lambda_electrostatics", "lambda_sterics"]
        values: list of float, like [1.0, 1.0, 1.0, 0.8, 0.5, 0.0], for different key(different lambda type), the len of list should be equal.
        Lambda states to calculate for each alchemical state
    simulation_lambdas : dict, have the similar structure with lambdas, but the values are the actual lambda values used for the simulation.
        Actual lambda values used for the simulation.
    simulation_lambdas_idxs : list
        Indices of lambda windows to simulate
    alchem_settings : utils.file_parser._SectionSettings
        Settings for the reconstructed non-bonded alchemical-alchemical interactions



    Notes
    -----
    When a PDBX file is provided, it enables group-based alchemical 
    transformations based on the specified atom sets and interaction groups.
    The class supports both classical FEP-ABFE and DC-ABFE alchemical transformations.
    """
    
    # pdbx is used to specify group vdw decouple and group charge annihilation
    def __init__(self, system_container, restraint_parm=None,
                 timestep=4*unit.femtoseconds, temperature=298*unit.kelvin,
                 pdbx=None,
                 current_group_nb=0, current_group_chg=0,
                 constant_pressure=None, plumed=None, read_state=None,
                 annihilate_electrostatics=True,
                 alchem_group='resname MOL',platform=None, platformProperties=None,
                 another_restraint_parm=None,
                 **kwargs):
        self.alchemical_state=None
        self.restraint_state=None
        self.timestep=timestep
        self.temperature = temperature
        self.context_cache=cache.ContextCache(platform=platform,platform_properties=platformProperties)
        self.original_system_container=system_container
        self.system_container=copy.deepcopy(system_container)
        self.restraint_parm=restraint_parm
        self.another_restraint_parm=another_restraint_parm
        self.read_state=read_state
        self.u_kln_array = None
        self.u_kln = None
        self.lambdas={}
        # added by rdliu
        self.simulation_lambdas={}
        self.simulation_lambdas_idxs=[]

        self.alchem_settings = kwargs['alchem_settings']  # Type: utils.file_parser._SectionSettings
        # try new setting class
        self.set_rbfe_exception = self.alchem_settings.set_rbfe_exception
        # added by ylzhong

        self.nsteps=0
        super().__init__(self.system_container)
        self.topology=self.original_system_container.topology

        self.current_group_nb=current_group_nb 

        self.pdbx=pdbx
        if pdbx is not None:
            self.pdbx = PdbxParser(pdbx)

            decouple_groups_nb = self.pdbx.get_group_nb_dict()
            initial_charge, target_charge = self.pdbx.get_charge_list(col=current_group_chg)
            # Use aligned atom id in pdbx.

            atomsets, interaction_groups = InteractionGroups(pdbx).run(self.alchem_settings.get_atomset_mode)
            self.alchem_settings.set_property('atomsets', atomsets)
            self.alchem_settings.set_property('interaction_groups', interaction_groups)
            # Use new method to get atomsets.

            print(decouple_groups_nb)
            print(initial_charge, target_charge)
        else:
            decouple_groups_nb=False
            initial_charge=False
            target_charge=False

        if plumed:
            self.useplumed(plumed)
        if constant_pressure:
            self.thermodynamic_state = states.ThermodynamicState(self.system_container.system, temperature, constant_pressure)
        else:
            self.thermodynamic_state = states.ThermodynamicState(self.system_container.system, temperature)

        self.alchemical_state = self.create_alchem_state(alchem_group, annihilate_electrostatics=annihilate_electrostatics,
                                                         set_initial_charge=initial_charge, set_target_charge=target_charge,
                                                         decouple_groups_nb=decouple_groups_nb, current_group_nb=current_group_nb,
                                                         alchem_settings=self.alchem_settings) #init self.alchemical_state
        # modded by ylzhong set_rbfe_exception
        # determine whether we need to add restraint
        if self.restraint_parm is not None:
            self.restraint_state = self.create_restrain_state(restraint_parm) #init self.restraint_state
            composable_states = [self.alchemical_state, self.restraint_state]
            if self.another_restraint_parm is not None:
                self.another_restraint_state = self.create_restrain_state(self.another_restraint_parm,True)
                composable_states.append(self.another_restraint_state)            
        else:
            composable_states = [self.alchemical_state]

        self.compound_state = states.CompoundThermodynamicState(
                         thermodynamic_state=self.thermodynamic_state, composable_states=composable_states)

        self.integrator= integrators.LangevinIntegrator(timestep=self.timestep,splitting="V R R R O R R R V")
        self.integrator = openmm.openmm.LangevinMiddleIntegrator(self.temperature, 1/unit.picoseconds, self.timestep)
        self.context, self.integrator = self.context_cache.get_context(self.compound_state, self.integrator)

        if not read_state:
            self.context.setPositions(self.system_container.positions)
        else:
            if isinstance(read_state, str):
                self.loadState(read_state) # change this to sampler_state
            else:
                self.set_state_xpv(read_state)
    
    def get_alchem_start_id(self, topology, alchem_group):
        """
        Get the starting atom index of the alchemical group. [Deprecated]
        
        Parameters
        ----------
        topology : mdtraj.Topology or openmm.Topology
            Molecular topology
        alchem_group : str
            DSL selection string for the alchemical group
            
        Returns
        -------
        int
            Index of the first atom in the alchemical group
        """
        if isinstance(topology, mdtraj.Topology):
            mdtraj_topology = topology
        else:
            mdtraj_topology = mdtraj.Topology.from_openmm(topology)
        # Determine indices of the atoms to restrain.
        alchem_atoms = mdtraj_topology.select(alchem_group).tolist()
        return alchem_atoms[0]

    def set_lambdas(self, lambdas):
        """
        Set lambda values for alchemical transformations.

        Parameters
        ----------
        lambdas : dict
            Dictionary containing lambda schedules for different components.
            All lambda components must have the same number of windows.

        Raises
        ------
        ValueError
            If the lengths of lambda schedules for different components are not equal.
        """
        self.lambdas=lambdas

        # assert the length of all lambdas to be the same
        len_prev_lambda = -1 # we do not have the length of previous lambda here
        for lambda_type in self.lambdas:
            len_current_lambda = len(self.lambdas[lambda_type])
            if len_prev_lambda > -1:
                if len_prev_lambda != len_current_lambda:
                    raise ValueError('length of lambdas should be the same')
            len_prev_lambda = len_current_lambda

    # added by rdliu
    def set_simulation_lambdas(self, lambdas):
        """
        Set lambda values for actual simulation windows and validate against MBAR computation lambdas.

        Parameters
        ----------
        lambdas : dict
            Dictionary containing lambda schedules for simulation windows.
            All lambda components must have the same number of windows.

        Raises
        ------
        ValueError
            If any of the following conditions are met:
            - The lengths of lambda schedules for different components are not equal
            - The lambda keys for simulation do not match those for MBAR computation
            - Any simulation lambda state is not present in the MBAR computation lambdas

        Notes
        -----
        This method:
        1. Validates that all lambda components have the same length
        2. Verifies that simulation lambda keys match MBAR computation lambda keys
        3. Ensures simulation lambdas are a subset of MBAR computation lambdas
        4. Creates indices mapping simulation lambdas to MBAR computation lambdas

        The simulation_lambdas_idxs attribute stores the indices of simulation 
        lambda states within the full set of MBAR computation lambda states.
        """
        self.simulation_lambdas = lambdas
        # assert the length of all lambdas to be the same
        len_prev_lambda = -1 # we do not have the length of previous lambda here
        for lambda_type in self.simulation_lambdas:
            len_current_lambda = len(self.simulation_lambdas[lambda_type])
            if len_prev_lambda > -1:
                if len_prev_lambda != len_current_lambda:
                    raise ValueError('length of lambdas should be the same')
            len_prev_lambda = len_current_lambda
        # get the keyname startwith "lambda_"
        compute_mbar_lambda_keys = [key for key in self.lambdas if key.startswith("lambda_")]
        simulation_lambda_keys = [key for key in self.simulation_lambdas if key.startswith("lambda_")]
        if compute_mbar_lambda_keys != simulation_lambda_keys:
            raise ValueError("Error! Your simulation_lambda_keys are not as same as the lambda_keys in the compute_mbar_lambda!")

        # judge if self.simulation_lambdas is the subset of self.lambdas
        compute_mbar_tuples = []
        simulation_tuples = []
        for i in range(len(self.lambdas[compute_mbar_lambda_keys[0]])):
            compute_tuple = tuple(self.lambdas[key][i] for key in compute_mbar_lambda_keys)
            compute_mbar_tuples.append(compute_tuple)
        for i in range(len(self.simulation_lambdas[compute_mbar_lambda_keys[0]])):
            simulation_tuple = tuple(self.simulation_lambdas[key][i] for key in compute_mbar_lambda_keys)
            simulation_tuples.append(simulation_tuple)   
        for item in simulation_tuples:
            # print(item, compute_mbar_tuples)
            if item not in compute_mbar_tuples:
                # print(item, compute_mbar_tuples)
                raise ValueError("Error! Your simulation lambda is not found in the mbar_lambda!")
            
        self.simulation_lambdas_idxs = [compute_mbar_tuples.index(item) for item in simulation_tuples]

    def change_state(self, state, nstate):
        """
        Update an alchemical state object with lambda values for a specified state.

        Parameters
        ----------
        state : AlchemicalState
            The alchemical state object to be modified
        nstate : int
            Index of the lambda state to apply

        Notes
        -----
        Updates all lambda parameters in the state object using values from 
        self.lambdas at the specified state index.
        """
        for lambda_type in self.lambdas:
            setattr(state, lambda_type, self.lambdas[lambda_type][nstate])

    def get_state_lambda(self, nstate):
        """
        Get lambda types and values for a specified state, with special handling
        for electrostatics and sterics, for recording the potential energy and the lambda state information
        in the output csv file. The electrostatics and sterics values are transformed as (1 - lambda) 
        rounded to 3 decimals for the gromacs style lambda representation.

        Parameters
        ----------
        nstate : int
            Index of the lambda state

        Returns
        -------
        lambda_types : list
            List of lambda parameter names
        lambda_values : list
            List of lambda values, with electrostatics and sterics values 
            transformed as (1 - lambda)

        Notes
        -----
        Special handling for certain lambda types:
        - lambda_electrostatics: returns (1 - lambda) rounded to 3 decimals
        - lambda_sterics: returns (1 - lambda) rounded to 3 decimals
        All other lambda types return their direct values.
        """
        lambda_types=[]
        lambda_values=[]
        for lambda_type in self.lambdas:
            lambda_value=self.lambdas[lambda_type][nstate]
            if lambda_type == 'lambda_electrostatics':
                lambda_value = np.around(1-lambda_value, decimals=3)
            elif lambda_type == 'lambda_sterics':
                lambda_value = np.around(1-lambda_value, decimals=3)
            lambda_types.append(lambda_type)
            lambda_values.append(lambda_value)
        return lambda_types, lambda_values

    def create_alchem_state(self, alchem_group, annihilate_electrostatics=True, set_initial_charge=False, set_target_charge=False, decouple_groups_nb=False, current_group_nb=0,
                            **kwargs):
        """
        Create an alchemical state for free energy calculations.

        Parameters
        ----------
        alchem_group : str
            MDTraj selection string defining the alchemical region
        annihilate_electrostatics : bool, optional
            Whether to annihilate electrostatics (default=True)
        set_initial_charge : bool, optional
            Whether to set initial charges for charge perturbation (default=False)
        set_target_charge : bool, optional
            Whether to set target charges for charge perturbation (default=False)  
        decouple_groups_nb : bool, optional
            Whether to use group-based nonbonded decoupling (default=False)
        current_group_nb : int, optional
            Current nonbonded group index for group-based decoupling (default=0)
        **kwargs
            Additional kwargs including:
            - alchem_settings: utils.file_parser._SectionSettings object for
            customizing the alchemical factory

        Returns
        -------
        alchemical_state : alchemy.AlchemicalState
            The created alchemical state object

        Notes
        -----
        Uses original_system_container instead of system_container for atom selection 
        to avoid potential exceptions.

        The method workflow:
        1. Selects alchemical atoms using MDTraj topology
        2. Creates AlchemicalRegion with specified parameters 
        3. Initializes AbsoluteAlchemicalFactory with optional settings
        4. Creates alchemical system with specified parameters
        5. Updates thermodynamic state with new alchemical system
        6. Generates and returns corresponding alchemical state
        """
        alchem_atoms = self.original_system_container.mdtraj_topology.select(alchem_group) # cannot use self.system_container, otherwise will raise an exception, which could be a bug
        alchemical_region = alchemy.AlchemicalRegion(alchemical_atoms=alchem_atoms,annihilate_electrostatics=annihilate_electrostatics)
        #print(alchemical_region.name)
        if 'alchem_settings' in kwargs.keys():
            alchem_settings = kwargs['alchem_settings']  # Type: utils.file_parser._SectionSettings
            factory = alchemy.AbsoluteAlchemicalFactory(alchem_settings=alchem_settings)
        else:
            factory = alchemy.AbsoluteAlchemicalFactory()

        alchemical_system = factory.create_alchemical_system(self.thermodynamic_state.system, alchemical_region,
                                                             set_initial_charge=set_initial_charge,
                                                             set_target_charge=set_target_charge,
                                                             decouple_groups_nb=decouple_groups_nb,
                                                             current_group_nb=current_group_nb,
                                                             )
        # modded by ylzhong 2024.10.29
        self.thermodynamic_state.system = alchemical_system
        alchemical_state = alchemy.AlchemicalState.from_system(alchemical_system)
        return alchemical_state

    def create_restrain_state(self, restraint, ifanother=False):
        """
        Create a restraint state using Boresch-style restraints. Supports define the second Boresch restraint for the double Boresch restraints.

        Parameters
        ----------
        restraint : object
            Object containing restraint parameters including:
            - kbond: bond force constant (kcal/mol/Å²)
            - ktheta1, ktheta2: angle force constants (kcal/mol/rad²)
            - kphi1, kphi2, kphi3: dihedral force constants (kcal/mol/rad²)
            - rec_atoms: list of receptor atoms to restrain
            - lig_atoms: list of ligand atoms to restrain
            - r: equilibrium bond distance (Å)
            - theta1, theta2: equilibrium angles (rad)
            - phi1, phi2, phi3: equilibrium dihedrals (rad)
        ifanother : bool, optional
            If True, creates a secondary restraint state using Boresch2 class
            (default=False)

        Returns
        -------
        restraint_state : RestraintState or RestraintState2
            The created restraint state object with lambda_restraints initialized to 0.0

        Notes
        -----
        Creates either a primary (Boresch) or secondary (Boresch2) restraint state
        with the following components:
        - One distance restraint
        - Two angle restraints
        - Three dihedral restraints

        All force constants are converted to OpenMM units:
        - Bonds: kcal/mol/Å² 
        - Angles/Dihedrals: kcal/mol/rad²

        The restraint is applied to the system's thermodynamic state before
        returning the restraint state object.
        """
        FORCE_CONSTANT_BOND=restraint.kbond*unit.kilocalories_per_mole/unit.angstrom**2
        K_theta1_FORCE_CONSTANT_ANGLE=restraint.ktheta1*unit.kilocalories_per_mole/unit.radians**2
        K_theta2_FORCE_CONSTANT_ANGLE=restraint.ktheta2*unit.kilocalories_per_mole/unit.radians**2
        K_phi1_FORCE_CONSTANT_DIHE=restraint.kphi1*unit.kilocalories_per_mole/unit.radians**2
        K_phi2_FORCE_CONSTANT_DIHE=restraint.kphi2*unit.kilocalories_per_mole/unit.radians**2
        K_phi3_FORCE_CONSTANT_DIHE=restraint.kphi3*unit.kilocalories_per_mole/unit.radians**2
        if not ifanother:
            restraint = Boresch(restrained_receptor_atoms=restraint.rec_atoms,
                                    restrained_ligand_atoms=restraint.lig_atoms,
                                    K_r=FORCE_CONSTANT_BOND, r_aA0=restraint.r*unit.angstrom,
                                    K_thetaA=K_theta1_FORCE_CONSTANT_ANGLE, K_thetaB=K_theta2_FORCE_CONSTANT_ANGLE,
                                    theta_A0=restraint.theta1*unit.radians, theta_B0=restraint.theta2*unit.radians,
                                    K_phiA=K_phi1_FORCE_CONSTANT_DIHE, K_phiB=K_phi2_FORCE_CONSTANT_DIHE, K_phiC=K_phi3_FORCE_CONSTANT_DIHE,
                                    phi_A0=restraint.phi1*unit.radians, phi_B0=restraint.phi2*unit.radians, phi_C0=restraint.phi3*unit.radians
            )
            restraint.restrain_state(self.thermodynamic_state)
            restraint_state = RestraintState(lambda_restraints=0.0)
        else:
            restraint = Boresch2(restrained_receptor_atoms=restraint.rec_atoms,
                                    restrained_ligand_atoms=restraint.lig_atoms,
                                    K_r=FORCE_CONSTANT_BOND, r_aA0=restraint.r*unit.angstrom,
                                    K_thetaA=K_theta1_FORCE_CONSTANT_ANGLE, K_thetaB=K_theta2_FORCE_CONSTANT_ANGLE,
                                    theta_A0=restraint.theta1*unit.radians, theta_B0=restraint.theta2*unit.radians,
                                    K_phiA=K_phi1_FORCE_CONSTANT_DIHE, K_phiB=K_phi2_FORCE_CONSTANT_DIHE, K_phiC=K_phi3_FORCE_CONSTANT_DIHE,
                                    phi_A0=restraint.phi1*unit.radians, phi_B0=restraint.phi2*unit.radians, phi_C0=restraint.phi3*unit.radians
            )
            restraint.restrain_state(self.thermodynamic_state)
            restraint_state = RestraintState2(lambda_restraints2=0.0)
        return restraint_state
        

    def run(self, nsteps, niterations, save_state=True, current_state=None, cal_all_ene=False, cal_adj_ene_num=5):
        """
        Run alchemical simulation with specified parameters.

        Parameters
        ----------
        nsteps : int
            Number of simulation steps per iteration
        niterations : int
            Total number of iterations to run for each state
        save_state : bool, optional
            Whether to save state files (default=True)
        current_state : int or None, optional
            Specific state to calculate. If None, calculates all states (default=None)
        cal_all_ene : bool, optional
            If True, calculates energy differences with respect to all states.
            If False, only calculates energy differences for adjacent states 
            (default=False)
        cal_adj_ene_num : int, optional
            Number of adjacent states to calculate energies for when cal_all_ene=False
            (default=5)

        Notes
        -----
        Workflow:
        1. Sets up simulation parameters and states based on lambda values
        2. For each state:
            - Creates progress bar
            - Checks for existing state files (state_s{state}.xml or state_g{group_nb}_s{state}.xml) to avoid redundant calculations
            - Determines which states to calculate energy differences for
            - Initializes pandas DataFrame for storing results
            - Runs iterations:
                - Updates system to current lambda state
                - Performs simulation steps
                - Calculates potential energies for relevant states
            - Saves state file and energy data to CSV

        Output files:
        - State files: state_[g{group_nb}_]s{state}.xml
        - Energy data: state_[g{group_nb}_]s{state}.csv

        Energy differences are calculated in units of kT, where:
        kT = NA * kB * temperature

        The CSV files contain energy differences organized with multi-index columns
        for lambda values and times.
        """
        # 
        # nsteps: run simulation for nsteps and then calculate deltaU, this is called one iteration
        # niteration: how many iterations will be calculated in total for each state
        # current_state: calculate a specific state, if not specified, all states will be calculated
        # cal_all_ene: If set to be True, then calculate the deltaU with respect to all states. 
        #              Otherwise only deltaU of adjacent states will be calculated.

        self.nsteps=nsteps
        kT = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * self.integrator.getTemperature()

        lambdas_per_state = self.alchem_settings.lambdas_per_state
        print(f"lambdas_per_state: {lambdas_per_state}")
        nstates = int(len(list(self.lambdas.values())[0]) / lambdas_per_state)
        # print(f"nstates: {nstates} {len(list(self.lambdas.values())[0])}")
        all_states = range(nstates)

        if lambdas_per_state != 1 and len(self.simulation_lambdas_idxs) != 0:
            self.simulation_lambdas_idxs = [int(i/lambdas_per_state)
                                            for i in self.simulation_lambdas_idxs if i % lambdas_per_state == 0]
        # modified by rdliu
        if not current_state:
            if len(self.simulation_lambdas_idxs) == 0:
                states = all_states
            else:
                states = self.simulation_lambdas_idxs       
        else: states = current_state

        state_skipped=False

        timestep = self.timestep.value_in_unit(unit=unit.femtoseconds)
        # print(f'states: {states}')
        for state in states:
            state_lambda_idx = int(state * lambdas_per_state)
            # progress bar
            try:
                from tqdm import tqdm
                # iterator = tqdm(range(niterations), ascii=True, desc="Alchemical state %3d" % (k))
                iterator = tqdm(range(niterations), ascii=True, desc="Alchemical state %3d" % (state),
                                unit='fs', unit_scale=self.nsteps*timestep)
            except:
                iterator = progressbar(range(niterations),           "Alchemical state %3d" % (state), 40)
            
            if save_state:
                if self.pdbx is not None:
                    prefix='state_g'+str(self.current_group_nb)+'_'
                else: prefix='state_'
                statefile=prefix+'s'+str(state)+'.xml'

            # skip state that has already finished
            if save_state and os.path.exists(statefile):
                last_state_file=statefile
                state_skipped=True
                print('State file %s already exists, skip this state', statefile)
                continue
            if state_skipped:
                self.loadState(last_state_file) 

            # calculate energy with respect to all lambdas or just adjacent lambdas
            if cal_all_ene:
                ene_states=all_states
            else:
                ene_states=get_adjacent_numbers(state_lambda_idx,
                                                int(nstates*lambdas_per_state) - 1, cal_adj_ene_num, lambdas_per_state)
            n_ene_states = len(ene_states)
            print(ene_states)

            #Initialize dU pandas dataframe
            simulation_lambda_types, simulation_lambda_values = self.get_state_lambda(state_lambda_idx)
            simulation_lambda_types.insert(0, 'times(ps)')
            # print(muti_idx_names)
            times_lambda_tuples = [(i*nsteps*timestep/1000,)+tuple(simulation_lambda_values)  for i in range(0, niterations)]
            muti_idx = pd.MultiIndex.from_tuples(times_lambda_tuples, names=simulation_lambda_types)
            zero_shape = np.zeros((niterations, n_ene_states))
            # columns_ = [tuple(self.get_state_lambda(l)[1]) for l in range(0,n_ene_states)]
            columns_ = [tuple(self.get_state_lambda(l)[1]) for l in ene_states] # fix bug for IndexError: iloc cannot enlarge its target object
            single_simulation_df = pd.DataFrame(zero_shape, columns=columns_)
            single_simulation_df.index = muti_idx
            # print(single_simulation_df.columns)
            # print(single_simulation_df.shape)
            # added by phologlucinol, to align the ene_states to the index of single_simulation_df
            index_gap = max(ene_states)+1-n_ene_states
            # start simulation for this state
            for iteration in iterator:
                self.change_state(self.compound_state, state_lambda_idx)
                self.compound_state.apply_to_context(self.context)
                self._simulate(endStep=self.currentStep+nsteps)
                # calcualte dU
                for l in ene_states:
                    if l > (all_states[-1] + 1) * lambdas_per_state - 1: continue
                    if l < all_states[0] * lambdas_per_state: continue
                    self.change_state(self.compound_state, l)
                    self.compound_state.apply_to_context(self.context)
                    single_simulation_df.iloc[iteration, l-index_gap] = self.context.getState(getEnergy=True).getPotentialEnergy() / kT

            # save state and dU files
            if save_state: self.saveState(statefile)

            if self.pdbx is not None:
                csv_prefix='state_g'+str(self.current_group_nb)+'_'
            else: csv_prefix='state_'
            csv_file=csv_prefix+'s'+str(state)+'.csv'
            single_simulation_df.to_csv(csv_file, sep="|")
            # np.save(dU_file, self.u_kln_array)

        #self.u_kln = self.format_u_kln(self.u_kln_array)

    @property
    def final_state(self):
        return self.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)


