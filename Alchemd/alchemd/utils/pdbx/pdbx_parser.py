from operator import ge
import numpy as np


class Atom():
    """
    Represents an atom in a molecular structure with its properties and coordinates.

    Parameters
    ----------
    atom_info : str
        String containing atom information in PDB-like format:
        "ATOM/HETATM atomid atomname resname resid x y z [charges] group_num"

    Attributes
    ----------
    atomid : int
        Atom identification number
    atomname : str
        Name of the atom
    resname : str
        Name of the residue containing the atom
    resid : int 
        Residue identification number
    coord : list of float
        XYZ coordinates [x, y, z]
    charge_info : list of float
        List of charge-related values
    group_num : str
        Group number identifier
    atomtype : str
        Atom type derived from atomname (excludes numbers)

    Methods
    -------
    update_chg_info(new_chg, col)
        Updates charge information at specified column index
    update_group_num(new_group_num)
        Updates the group number
    write_atm_line()
        Generates formatted PDB-style atom line
    
    Raises
    ------
    ValueError
        If atom_info string format is incorrect
    """
    def __init__(self, atom_info):
        info_list = atom_info.strip().split()
        len_info = len(info_list)
        if len_info >= 8 and (
                info_list[0] == 'ATOM' or info_list[0] == 'HETATM'):

            self.atomid = int(info_list[1])
            self.atomname = info_list[2]
            self.resname = info_list[3]
            self.resid = int(info_list[4])
            self.coord = [float(crd) for crd in info_list[5:8]]
            self.charge_info = [float(chg) for chg in info_list[8:-1]]
            self.group_num = info_list[-1]
        else:
            raise ValueError("Atom information not correct: %s" % (atom_info))

    @property
    def atomtype(self):
        return ''.join([i for i in self.atomname if not i.isdigit()])

    def update_chg_info(self, new_chg, col):
        """
        Updates charge information at specified column.

        Parameters
        ----------
        new_chg : float
            New charge value
        col : int
            Column index (1-based) in charge_info to update

        Raises
        ------
        AssertionError
            If new_chg is not float type
        """
        assert isinstance(new_chg, (float, np.float64)
                          ), "new_chg must be float"
        self.charge_info[col - 1] = new_chg

    def update_group_num(self, new_group_num,):
        """
        Updates the group number identifier.

        Parameters
        ----------
        new_group_num : int
            New group number value

        Raises
        ------
        AssertionError
            If new_group_num is not int type
        """
        assert isinstance(new_group_num, int), "new_group_num must be int"
        self.group_num = new_group_num

    def write_atm_line(self):
        """
        Generates a formatted PDB-style atom line.

        Returns
        -------
        str
            Formatted atom line containing all atom information in PDB format
        """
        empty_altLoc = ' '
        empty_chainID = ' '
        empty_iCode = ' '
        formatted_charge_ = ["{: >10.6f}".format(
            chg) for chg in self.charge_info]
        charge_info_str = " ".join(formatted_charge_)
        nb_group_idx = self.group_num
        # print(charge_info_str)
        l_str = '{: <6s}{: >5d} {: ^4s}{}{:3s} {}{:4d}{}   {:8.3f}{:8.3f}{:8.3f} {} {}\n'.format(
            'HETATM',
            self.atomid,
            self.atomname,
            empty_altLoc,
            self.resname,
            empty_chainID,
            self.resid,
            empty_iCode,
            self.coord[0],
            self.coord[1],
            self.coord[2],
            charge_info_str,
            nb_group_idx,
        )
        return l_str

    def __repr__(self):
        return "Atom('" + ':%-2s@%-2s' % (self.resname, self.atomname) + "')"

    def __str__(self):
        return self.atomname


class PdbxParser():
    """
    Parser for PDBx format files that handles atom properties and group-based transformations.

    Parameters
    ----------
    pdbx_file : str
        Path to input PDBx format file

    Attributes
    ----------
    pdbx : list
        Raw lines from PDBx file
    atoms_list : list
        List of Atom objects parsed from PDBx file

    Examples
    --------
    Basic usage:
    >>> pdbx_obj = PdbxParser("input.pdbx")
    >>> grp_dict = {"group_1": [1,2,3], "group_2": [4,5,6]} 
    
    Decharge all atoms, VDW by groups:
    >>> pdbx_obj.dechg_all_devdw_bygroups(grp_dict)
    >>> pdbx_obj.writePDBX("output.pdbx")

    Group-based charge and VDW annihilation:
    >>> pdbx_obj.annihilate_grps(grp_dict, 'balanced')
    >>> pdbx_obj.writePDBX("output.pdbx")

    Notes
    -----
    The class provides different strategies for modifying atom properties:
    - Full decharging with grouped VDW via dechg_all_devdw_bygroups()
    - Group-based charge and VDW annihilation via annihilate_grps()
    
    Modified PDBx files can be written using writePDBX() method.
    """
    def __init__(self, pdbx_file):
        self.pdbx = self.loadPDBX(pdbx_file)
        self.atoms_list = []
        for line in self.pdbx:
            try:
                atom = Atom(line)
                self.atoms_list.append(atom)
            except ValueError:
                pass

    def loadPDBX(self, file):
        """
        Loads and reads PDBx format file.

        Parameters
        ----------
        file : str or file object
            Path to PDBx file or file object

        Returns
        -------
        list
            Lines read from PDBx file
        """
        if isinstance(file, str):
            with open(file, 'r') as f:
                pdbx = f.readlines()
        else:
            pdbx = file.readlines()
        return pdbx

    def writePDBX(self, file_name):
        """
        Writes atoms(Atom class) to a PDBx format file.

        Parameters
        ----------
        file_name : str
            Output file path
        """
        with open(file_name, 'w', encoding='utf-8') as f:
            for atom in self.atoms_list:
                f.write(atom.write_atm_line())

    # def get_group_nb_dict(self,start_id):
    #     group_nb={}
    #     atom_id=0
    #     for atom in self.atoms_list:
    #         group_name='group'+str(atom.group_num)
    #         if group_name in group_nb:
    #             group_nb[group_name].add(atom_id+start_id)
    #         else:
    #             group_nb[group_name]={atom_id+start_id} # create this set
    #         atom_id += 1
    #     return group_nb
    def get_group_nb_dict(self):
        """
        Creates dictionary mapping group names to atom IDs.

        Returns
        -------
        dict
            Dictionary with keys 'group<number>' and values as sets of atom IDs
        """
        group_nb = {}
        for atom in self.atoms_list:
            group_name = 'group' + str(atom.group_num)
            if group_name in group_nb:
                # since the atomid in openmm is started from 0
                group_nb[group_name].add(atom.atomid-1)
            else:
                group_nb[group_name] = {atom.atomid-1}  # create this set
        return group_nb

    def get_charge_list(self, col=1):
        """
        Extracts initial and target charges from specified columns.

        Parameters
        ----------
        col : int, default=1
            Column index (1-based) for charge information.
            If 0, uses same column for initial and target charges.

        Returns
        -------
        tuple of dict
            Two dictionaries containing:
            - Initial charges mapped to atom IDs (0-based)
            - Target charges mapped to atom IDs (0-based)
            
        Notes
        -----
        If target charge column doesn't exist, target charges default to 0.0
        """
        initial_charge_dict = {}
        target_charge_dict = {}
        # atom_id=0
        if col == 0:
            initial_col = col
            target_col = col
        else:
            initial_col = col - 1
            target_col = col
        for atom in self.atoms_list:
            # initial_charge_dict[atom_id+start_id]=atom.charge_info[initial_col]
            initial_charge_dict[atom.atomid -
                                1] = atom.charge_info[initial_col]
            try:
                # target_charge_dict[atom_id+start_id]=atom.charge_info[target_col]
                target_charge_dict[atom.atomid -
                                   1] = atom.charge_info[target_col]
            except BaseException:
                # target_charge_dict[atom_id+start_id]=0.0
                target_charge_dict[atom.atomid - 1] = 0.0
            # atom_id += 1
        return initial_charge_dict, target_charge_dict

    def get_net_chg_of_atms_by_atmid_lst(self, atomid_list, col):
        """
        Calculates net charge for specified atoms at given charge column.

        Parameters
        ----------
        atomid_list : list
            List of atom IDs from PDBx file
        col : int
            Column index (1-based) in atom.charge_info

        Returns
        -------
        float
            Sum of charges for specified atoms
        """
        net_chg = 0
        for atomid in atomid_list:
            for atom in self.atoms_list:
                if atom.atomid == atomid:
                    net_chg += atom.charge_info[col - 1]
        return net_chg

    def find_positive_chg_atms(self, atomid_list, which_chg_col):
        """
        Finds atoms with positive charges from specified atom IDs.

        Parameters
        ----------
        atomid_list : list
            List of atom IDs from PDBx file
        which_chg_col : int
            Column index (1-based) in atom.charge_info

        Returns
        -------
        list
            List of Atom objects having positive charge in specified column
        """
        posi_chg_atms = []
        for atm in self.atoms_list:
            if atm.charge_info[which_chg_col -
                               1] > 0 and atm.atomid in atomid_list:
                posi_chg_atms.append(atm)
        return posi_chg_atms

    def find_negative_chg_atms(self, atomid_list, which_chg_col):
        """
        Finds atoms with negative charges from specified atom IDs.

        Parameters
        ----------
        atomid_list : list
            List of atom IDs from PDBx file
        which_chg_col : int
            Column index (1-based) in atom.charge_info

        Returns
        -------
        list
            List of Atom objects having negative charge in specified column
        """
        nega_chg_atms = []
        for atm in self.atoms_list:
            if atm.charge_info[which_chg_col -
                               1] < 0 and atm.atomid in atomid_list:
                nega_chg_atms.append(atm)
        return nega_chg_atms

    def update_atom_prop(self, atomid, prop_name, new_value, **kwargs):
        """
        Updates property of specified atom.

        Parameters
        ----------
        atomid : int
            Atom ID from PDBx file
        prop_name : str
            Name of property to update
        new_value : any
            New value for the property
        **kwargs : dict
            Additional keyword arguments:
            - which_chg_col : int
                Column index (1-based) in atom.charge_info to update

        Notes
        -----
        If which_chg_col is provided, updates specific charge in charge_info list.
        Otherwise directly sets the specified property to new_value.
        """
        for atm in self.atoms_list:
            if atm.atomid == atomid:
                if 'which_chg_col' in kwargs.keys():
                    old_chg_info = getattr(atm, prop_name)
                    old_chg_info[kwargs['which_chg_col'] - 1] = new_value

                    setattr(atm, prop_name, old_chg_info)
                else:
                    setattr(atm, prop_name, new_value)

    def get_atom_prop(self, atomid, prop_name):
        """
        Gets property value for specified atom.

        Parameters
        ----------
        atomid : int
            Atom ID from PDBx file
        prop_name : str
            Name of property to retrieve

        Returns
        -------
        any
            Value of requested property for specified atom
        """
        for atm in self.atoms_list:
            if atm.atomid == atomid:
                prop = getattr(atm, prop_name)
        return prop

    def dup_last_chg_info_col(self, ):
        """
        Duplicates last charge information column for all atoms.
        Appends copy of last charge value to charge_info list of each atom.
        """
        for atm in self.atoms_list:
            atm.charge_info.append(atm.charge_info[-1])

    def annihilate_grps(self, grp_dict, dechg_strategy='balanced'):
        """
        Applies group-based charge annihilation.

        Parameters
        ----------
        grp_dict : dict
            Dictionary mapping group names to atom ID lists
            Format: {'group_1': [atomid_1, atomid_2, ...], ...}
        dechg_strategy : str, default='balanced'
            Strategy for decharging atoms. Currently supports:
            - 'balanced': Uses balanced decharging approach

        Notes
        -----
        Assigns sequential group numbers starting from 1 to atoms in each group.
        """
        start_grp_id = 1
        for key, value in grp_dict.items():
            if dechg_strategy == 'balanced':
                self.balanced_dechg_one_grp(value, start_grp_id)
            for atom in self.atoms_list:
                if atom.atomid in value:
                    atom.group_num = start_grp_id
            start_grp_id += 1

    def dechg_all_devdw_bygroups(self, grp_dict):
        """
        Simultaneously decharges all atoms and decouples vdW by groups.

        Parameters
        ----------
        grp_dict : dict
            Dictionary mapping group names to atom ID lists
            Format: {'group_1': [atomid_1, atomid_2, ...], ...}

        Notes
        -----
        - Sets charges to [0.0, 0.0] for all atoms
        - Assigns sequential group numbers for vdW decoupling
        """
        start_grp_id = 1
        for key, value in grp_dict.items():
            for atom in self.atoms_list:
                if atom.atomid in value:
                    atom.group_num = start_grp_id
                    atom.charge_info = [0.0, 0.0]
            start_grp_id += 1

    def dechg_bygroup_devdw_bygroups(self, grp_dict):
        """
        Decharges atoms group-by-group and decouples vdW by groups.

        Parameters
        ----------
        grp_dict : dict
            Dictionary mapping group names to atom ID lists
            Format: {'group_1': [atomid_1, atomid_2, ...], ...}

        Notes
        -----
        - Assigns sequential group numbers
        - Duplicates charge columns if needed
        - Sets charges to 0 for group 2 atoms in new charge column
        """
        start_grp_id = 1
        for key, value in grp_dict.items():
            for atom in self.atoms_list:
                if atom.atomid in value:
                    atom.group_num = start_grp_id
                    example_atom_chg_info = self.get_atom_prop(
                        self.atoms_list[0].atomid, 'charge_info')
                    if len(example_atom_chg_info) == start_grp_id:
                        self.dup_last_chg_info_col()
                    if atom.group_num == 2:
                        self.update_atom_prop(
                            atom.atomid, 'charge_info', 0.0000, which_chg_col=start_grp_id + 1)
            start_grp_id += 1

    def balanced_dechg_one_grp(self, atomid_list, bef_chg_info_col, ):
        """
        Performs balanced decharging of one group while maintaining system net charge.

        Parameters
        ----------
        atomid_list : list
            List of atom IDs from PDBx file to be decharged
        bef_chg_info_col : int 
            Column index (1-based) in atom.charge_info containing charges to be changed
            The target column will be bef_chg_info_col + 1

        Notes
        -----
        Algorithm:
        1. Duplicates last charge column if needed
        2. Sets charges to 0 for specified atoms in target column
        3. Redistributes the net charge change to remaining atoms:
            - If net charge was negative: distributed to positive atoms
            - If net charge was positive: distributed to negative atoms
        4. Distribution is weighted by original charges
        
        Charges are rounded to 6 decimal places.
        
        Prints:
        - Initial net charge of group before changes
        - Final net charge of system after balanced redistribution

        Example
        -------
        If group had -1.0 net charge:
        - Group atoms set to 0
        - -1.0 charge distributed among remaining positive atoms
        proportionally to their original charges to maintain neutrality
        """
        aft_chg_info_col = bef_chg_info_col + 1
        example_atom_chg_info = self.get_atom_prop(
            self.atoms_list[0].atomid, 'charge_info')
        if len(example_atom_chg_info) == bef_chg_info_col:
            self.dup_last_chg_info_col()
        before_change_grp_netchg = np.around(
            self.get_net_chg_of_atms_by_atmid_lst(
                atomid_list, bef_chg_info_col), decimals=6)
        print(f'before_change_grp_netchg: {before_change_grp_netchg}')
        # set the decharge group's aft_chg_info
        for atomid in atomid_list:
            self.update_atom_prop(
                atomid,
                'charge_info',
                0.0000,
                which_chg_col=aft_chg_info_col)
        # get the remaining non-zero-charge atomid_list
        remain_non_zero_atomid_list = []
        for atm in self.atoms_list:
            # first check if the atom is in the atomid_list that we want to
            # annihilate
            if atm.atomid not in atomid_list:
                # second check if the atom's bef_chg_info is non-zero
                if atm.charge_info[bef_chg_info_col] != 0.0:
                    remain_non_zero_atomid_list.append(atm.atomid)
        # now assign the before_change_grp_netchg to the remaining non-zero-charge atom
        # if the remain_non_zero_atomid_list is empty list, then skip the
        # assignment
        if len(remain_non_zero_atomid_list) == 0:
            pass
        else:
            if before_change_grp_netchg < 0:
                positive_chg_atms = self.find_positive_chg_atms(
                    remain_non_zero_atomid_list, bef_chg_info_col)
                positive_chg_atmsid_list = [
                    atm.atomid for atm in positive_chg_atms]
                net_posi_chg = np.around(
                    self.get_net_chg_of_atms_by_atmid_lst(
                        positive_chg_atmsid_list,
                        bef_chg_info_col),
                    decimals=6)
                for atomid in positive_chg_atmsid_list:
                    # weighted assignment
                    old_chg = self.get_atom_prop(atomid, 'charge_info')[
                        bef_chg_info_col - 1]
                    new_chg = np.around(
                        old_chg +
                        before_change_grp_netchg *
                        old_chg /
                        net_posi_chg,
                        decimals=6)
                    self.update_atom_prop(
                        atomid, 'charge_info', new_chg, which_chg_col=aft_chg_info_col)
            else:
                negative_chg_atms = self.find_negative_chg_atms(
                    remain_non_zero_atomid_list, bef_chg_info_col)
                negative_chg_atmsid_list = [
                    atm.atomid for atm in negative_chg_atms]
                net_nega_chg = np.around(
                    self.get_net_chg_of_atms_by_atmid_lst(
                        negative_chg_atmsid_list,
                        bef_chg_info_col),
                    decimals=6)
                for atomid in negative_chg_atmsid_list:
                    old_chg = np.around(
                        self.get_atom_prop(
                            atomid, 'charge_info')[
                            bef_chg_info_col - 1], decimals=6)
                    new_chg = np.around(
                        old_chg +
                        before_change_grp_netchg *
                        old_chg /
                        net_nega_chg,
                        decimals=6)
                    self.update_atom_prop(
                        atomid, 'charge_info', new_chg, which_chg_col=aft_chg_info_col)
        all_atom_ids_list = [atom.atomid for atom in self.atoms_list]
        after_dechg_sys_netchg = self.get_net_chg_of_atms_by_atmid_lst(
            all_atom_ids_list, aft_chg_info_col)
        print(
            f'After balanced charge assignment, the net charge of the ligand is {after_dechg_sys_netchg}')


# added by zdli
# To get interaction groups for alchemically_modify_NonbondedForce
class InteractionGroups:
    """
    Class for determining interaction groups for alchemical nonbonded force modifications.

    Parameters
    ----------
    pdbx_file : str
        Path to PDBx format input file

    Notes
    -----
    Supports three interaction group modes:
    - normal_by_group: Single alchemical region interactions
    - cs_fep: combined-structure FEP (Reference: Z. Li, M.-Y. Jiang, R. Liu, Q. Wang, Q. Zhou, Y.-Y. Huang, Y. Wu, C.-G. Zhan and H.-B. Luo, Discovery of Highly Potent Phosphodiesterase-1 Inhibitors by a Combined-Structure Free Energy Perturbation Approach, Acta Pharmaceutica Sinica B, 2024, S2211383524002521, DOI: 10.1016/j.apsb.2024.06.021.)
    - pseudocore: Pseudocore-based transformations
    """
    def __init__(self, pdbx_file):

        self.pdbx = PdbxParser(pdbx_file)

    def _get_atomdict_with_resname(self, resname: str) -> dict:
        """
        Gets dictionary mapping atom IDs to group numbers for a residue.

        Parameters
        ----------
        resname : str
            Residue name to filter atoms

        Returns
        -------
        dict
            Dictionary mapping atom IDs to group numbers {atom.atomid: atom.group_num}
        """
        return_atom_dict = {}
        for atom in self.pdbx.atoms_list:
            if atom.resname == resname:
                return_atom_dict[atom.atomid] = int(atom.group_num)
        return return_atom_dict

    def normal_by_group(self) -> (tuple, set):
        """
        Creates interaction groups for single molecule transformations.

        Returns
        -------
        tuple
            Single frozenset containing all molecule atom IDs
        set
            Set containing single tuple of (mol_atomset, mol_atomset) for self-interactions
        """
        mol_atomset = set()
        for atom in self.pdbx.atoms_list:
            if atom.resname == 'MOL':
                mol_atomset.add(atom.atomid-1)
        mol_atomset = frozenset(mol_atomset)
        atomsets = (mol_atomset)
        interaction_groups = {(mol_atomset, mol_atomset)}
        return atomsets, interaction_groups

    def common_structure_fep(self) -> (tuple, set):
        """
        Creates interaction groups for common structure FEP.

        Returns
        -------
        tuple
            (unit_A_atomset, unit_B_atomset, common_atomset) where:
            - unit_A_atomset: Unique atoms in state A (LAM residue)
            - unit_B_atomset: Unique atoms in state B (LBM residue) 
            - common_atomset: Shared atoms (COM residue)
        set
            Set of interaction tuples defining allowed interactions between atom sets

        Raises
        ------
        ValueError
            If atom assignments invalid:
            - Overlapping atom sets
            - Multiple discontinuous groups in both end states
            - Multiple groups in common atoms
        """

        # LAM means unique atoms of A,LBM means unique atoms of B; COM means
        # shared atoms of A and B.
        unit_A_atomdict = self._get_atomdict_with_resname('LAM')
        unit_B_atomdict = self._get_atomdict_with_resname('LBM')
        common_atomdict = self._get_atomdict_with_resname('COM')

        unit_A_atomset = frozenset(i-1 for i in unit_A_atomdict.keys())
        unit_B_atomset = frozenset(i-1 for i in unit_B_atomdict.keys())
        common_atomset = frozenset(i-1 for i in common_atomdict.keys())

        if unit_A_atomset.intersection(unit_B_atomset):
            raise ValueError('Existing atoms unassigned to common_atomset!')
        elif unit_A_atomset.intersection(common_atomset):
            raise ValueError(
                'Please Check if only unique atoms or all unique atoms in LAM !')
        elif unit_B_atomset.intersection(common_atomset):
            raise ValueError(
                'Please Check if only unique atoms or all unique atoms in LBM !')

        # check if atomdict_grp_num meets requirement.
        # Only one of the ligand_atomsets can have more than one group.

        tobe_check_atomdict = [unit_A_atomdict, unit_B_atomdict]
        count = 0
        for atomdict in tobe_check_atomdict:
            if atomdict:
                atomdict_grp_num = set(atomdict.values())
                if len(atomdict_grp_num) > 1:
                    # check if the groups are continuous.
                    if len(atomdict_grp_num) != max(atomdict_grp_num) - min(atomdict_grp_num) + 1:
                        raise ValueError('Uncontinuous groups!')

                    count = count + 1
        if count == 2:
            raise ValueError(
                'Both unit_A and unit_B have more than one group! ')

        if len(set(common_atomdict.values())) > 1:
            raise ValueError('Common atoms must be in the same group!')

        atomsets = (unit_A_atomset, unit_B_atomset, common_atomset)
        interaction_groups = {(unit_A_atomset, unit_A_atomset),
                              (unit_B_atomset, unit_B_atomset),
                              (common_atomset, common_atomset),
                              (unit_A_atomset, common_atomset),
                              (unit_B_atomset, common_atomset)}

        print('Interaction groups:', interaction_groups)
        print('atomsets:', atomsets)

        return atomsets, interaction_groups

    def pseudocore(self) -> (tuple, set):
        """
        Creates interaction groups for pseudocore transformations.

        Returns
        -------
        tuple
            (ligand_A_atomset, ligand_B_atomset) containing all atoms in each end state
        set
            Set of interaction tuples for self-interactions of each end state

        Raises
        ------
        ValueError
            If atom assignments invalid:
            - No unique atoms in state B
            - Multiple discontinuous groups in both end states
        """

        # LAM means unique atoms of A,LBM means unique atoms of B; COM means
        # shared atoms of A and B.
        unit_A_atomdict = self._get_atomdict_with_resname('LAM')
        unit_B_atomdict = self._get_atomdict_with_resname('LBM')
        common_atomdict = self._get_atomdict_with_resname('COM')

        ligand_A_atomdict = {
            **unit_A_atomdict,
            **common_atomdict}
        ligand_B_atomdict = {
            **unit_B_atomdict,
            **common_atomdict}

        ligand_A_atomset = frozenset(i-1 for i in ligand_A_atomdict.keys())
        ligand_B_atomset = frozenset(i-1 for i in ligand_B_atomdict.keys())

        if len(ligand_A_atomset.union(ligand_B_atomset) - ligand_A_atomset) == 0:
            raise ValueError('No unique atoms in ligand_B !')

        # check if atomdict_grp_num meets requirement.
        # Only one of the ligand_atomsets can have more than one group.

        tobe_check_atomdict = [ligand_A_atomdict, ligand_B_atomdict]
        count = 0
        for atomdict in tobe_check_atomdict:

            if atomdict:
                atomdict_grp_num = set(atomdict.values())

                if len(atomdict_grp_num) > 2:
                    # check if the groups are continuous.
                    if len(atomdict_grp_num) != max(atomdict_grp_num) - min(atomdict_grp_num) + 1:
                        raise ValueError('Uncontinuous groups!')

                    count = count + 1
        if count == 2:
            raise ValueError(
                'Both ligand_A and ligand_B have more than two group! ')

        atomsets = (ligand_A_atomset, ligand_B_atomset, ())
        interaction_groups = {(ligand_A_atomset, ligand_A_atomset),
                              (ligand_B_atomset, ligand_B_atomset)}

        return atomsets, interaction_groups


    def run(self, mode):
        """
        Executes interaction group creation for specified mode.

        Parameters
        ----------
        mode : str
            Interaction mode:
            - 'normal_by_group': Single molecule
            - 'cs_fep': combined-structure FEP  
            - 'pseudocore': Pseudocore transformation

        Returns
        -------
        tuple
            Mode-specific atom sets
        set
            Mode-specific interaction groups
        """
        if mode == 'normal_by_group':
            atomsets, interaction_groups = self.normal_by_group()
        if mode == 'cs_fep':
            atomsets, interaction_groups = self.common_structure_fep()
        if mode == 'pseudocore':
            atomsets, interaction_groups = self.pseudocore()
        return atomsets, interaction_groups




