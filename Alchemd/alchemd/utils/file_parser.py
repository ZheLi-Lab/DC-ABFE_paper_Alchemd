import types
import warnings
'''
###Input file format:
[md_control]
temperature = 300
time_step = 0.01
num_steps = 10000

[md_system]
mass=10

[pme]
use_pme = True
pme_cutoff =12.0

###To use this class:
# Create an instance of the InputParser class
input_parser = InputParser('input.txt')

# Get the parsed data for the md_control section
md_control = input_parser.get_md_control()

# Get the parsed data for the md_system section
md_system = input_parser.get_md_system()

# Get the parsed data for the pme section
pme = input_parser.get_pme()

# Use the parsed data in your molecular dynamics program
run_simulation(md_control, md_system, pme)
'''


class _SectionSetting:
    def __init__(self, name):
        self.name = name
        self.__prop_names = set()

    def append_to_file(self, file):
        with open(file, 'a') as f:
            f.write(f'{self.__str__()}')

    def set_property(self, key: str, value):
        """
        _SectionSetting.set_property(key, value), the prop is read-only.
        If you want to change the property, use set_property().
        You can pass the property by using Instance.key.
        :param key: the key of the property
        :param value: the value of the property
        :return: None
        """
        self.__prop_names.add(key)
        if f"_{key}" in self.__dict__.keys():
            warnings.warn(f"SectionSetting {self.name} already has an attribute <{key}>")
        setattr(self, f"_{key}", value)
        setattr(self.__class__, key, property(lambda _self, _key=key: getattr(_self, f"_{_key}")))

    def get_property(self, key: str):
        return getattr(self, f"_{key}")

    def has_property(self, key: str):
        return key in self.__prop_names

    def __str__(self):
        return_block = f'[{self.name}]\n'
        for key in self.__prop_names:
            return_block += f'{key} = {self.get_property(key)}\n'
        return return_block

    def apply_with_dict(self, props: dict):
        for key, value in props.items():
            self.set_property(key, value)


class AllSetting:
    def __init__(self):
        self.section_names = []

    def add_section(self, name: str, section_setting: _SectionSetting):
        if isinstance(section_setting, _SectionSetting):
            self.section_names.append(name)
            setattr(self, f"_{name}", section_setting)
            setattr(self.__class__, name, property(lambda _self, _name=name: getattr(_self, f"_{_name}")))
        else:
            raise TypeError(f"Settings of {name} is not a SectionSettings")

    def get_section(self, name: str) -> _SectionSetting:
        return getattr(self, f"_{name}")

    def has_section(self, name: str) -> bool:
        return name in self.section_names

    def __str__(self):
        return_block = ''
        for name in self.section_names:
            return_block += f"{self.get_section(name)}"
        return return_block

    def apply_with_dict(self, props: dict):
        """
        the dict should contain section setting dict.
        Format: {'section_name': {'section_setting': section_setting,},}
        :param props:
        :return:
        """
        for name, section_setting in props.items():
            getattr(self, f"_{name}").apply_with_dict(section_setting)

    def autocomplete_with_default(self, props: dict):
        for section_name, section_settings in props.items():
            if self.has_section(section_name):
                for setting_name, setting_value in section_settings.items():
                    if not self.get_section(section_name).has_property(setting_name):
                        self.get_section(section_name).set_property(setting_name, setting_value)
            else:
                tmp_section = _SectionSetting(section_name)
                for setting_name, setting_value in section_settings.items():
                    tmp_section.set_property(setting_name, setting_value)
                self.add_section(section_name, tmp_section)



class InputParser:
    def __init__(self, filename):
        # Define the expected keys for each section and their default values
        self.def_sections()

        # Auto add props to self.sections. Add by ylzhong.
        self.dynamic_section_mod(filename)

        # Initialize variables to hold the parsed data
        self.data = {}
        for section, options in self.sections.items():
            self.data[section] = {}
            for option, default_value in options.items():
                self.data[section][option] = default_value

        self.current_section = None
        tmp_section_setting = _SectionSetting('emp')

        self.all_setting = AllSetting()  # Add by ylzhong

        # Parse the input file
        with open(filename, 'r') as f:
            for line in f:
                # Remove inline comments
                line = line.split('#')[0].strip()
                if not line:  # Skip empty lines and comments
                    continue
                if line.startswith('[') and line.endswith(']'):
                    if self.current_section is not None:
                        self.all_setting.add_section(self.current_section, tmp_section_setting)
                    self.current_section = line[1:-1].lower()
                    tmp_section_setting = _SectionSetting(self.current_section)
                elif self.current_section is not None:
                    for item in line.split(','):
                        key, value = item.split('=')
                        key = key.strip().lower()
                        value = value.strip()
                        if key in self.sections[self.current_section]:
                            if value.isdigit():
                                formatted_value = int(value)
                            elif value.replace('.', '').isdigit():
                                formatted_value = float(value)
                            elif value.lower() == 'true':
                                formatted_value = True
                            elif value.lower() == 'false':
                                formatted_value = False
                            elif value.lower() == 'none':
                                formatted_value = None
                            else:
                                formatted_value = value.strip('"').strip("'")

                            self.data[self.current_section][key] = formatted_value
                            tmp_section_setting.set_property(key, formatted_value)

        self.all_setting.add_section(self.current_section, tmp_section_setting)
        # print(self.all_setting)
        self.all_setting.autocomplete_with_default(self.sections)

        # Define the section-related methods dynamically
        for section in self.sections:
            def get_section_data(self, section=section):
                return self.data.get(section, {})
            setattr(self, f'get_{section}', types.MethodType(get_section_data, self))

    def get_all_settings(self):
        return self.all_setting

    def dynamic_section_mod(self, input_file):
        """
        Dynamically load the section props' names from the input file.
        Adding new props to the setting is no longer required to change self.sections manually.
        :return:
        """
        tmp_now_section_name = None
        tmp_now_section_setting_dict = {}
        tmp_all_section_setting_dict = {}
        with open(input_file, 'r') as f:
            for line in f:
                # Remove inline comments
                line = line.split('#')[0].strip()
                if not line:  # Skip empty lines and comments
                    continue
                if line.startswith('[') and line.endswith(']'):
                    if tmp_now_section_name is not None:
                        tmp_all_section_setting_dict[tmp_now_section_name] = tmp_now_section_setting_dict
                    tmp_now_section_name = line[1:-1].lower()
                    tmp_now_section_setting_dict = {}
                elif tmp_now_section_name is not None:
                    for item in line.split(','):
                        key, value = item.split('=')
                        key = key.strip().lower()
                        value = value.strip()

                        if value.isdigit():
                            formatted_value = int(value)
                        elif value.replace('.', '').isdigit():
                            formatted_value = float(value)
                        elif value.lower() == 'true':
                            formatted_value = True
                        elif value.lower() == 'false':
                            formatted_value = False
                        elif value.lower() == 'none':
                            formatted_value = None
                        else:
                            formatted_value = value.strip('"').strip("'")

                        tmp_now_section_setting_dict[key] = formatted_value

            tmp_all_section_setting_dict[tmp_now_section_name] = tmp_now_section_setting_dict
            # print(tmp_all_section_setting_dict)
            for key, section_dict in tmp_all_section_setting_dict.items():
                if key not in self.sections:
                    self.sections[key] = section_dict
                    continue

                for prop, value in section_dict.items():
                    if prop not in self.sections[key].items():
                        self.sections[key][prop] = value

    def def_sections(self):
        # Define the expected keys for each section and their default values
        self.sections = {
            # 'job_type': {
            #     'jobtype': 'normal' # if job type is normal[default], then use following options to control a single job.
            #     },            # TODO: FEP-Cascade, customized job type that easily perform FEP calculations.
            'normalmd': {
                'normalmd' : False, # whether perform normal MD
                'temperature': 298.0, 
                'timestep': 2,  #unit in fs, could be 4 fs, but it may cause unstablity of MD
                'nsteps': 1000,
                'save_traj': False
                },
            'alchemical': {
                'alchemical': False, #whether to run alchemical MD
                'mode': 'normal',
                'lambdas_json' : 'lambdas.json',
                'lambdas_group' : 'lambda_com_21normal',
                'simulation_lambdas_name': None,
                'set_rbfe_exception': False,
                'get_atomset_mode': 'cs_fep',  # normal_by_group, cs_fep, pseudocore
                'softcore': 'ssc2',  # gapsys, ssc2 or old
                'lambdas_per_state': 1,  # how many lambdas in lambda_json count as a single state.
                'temperature': 298.0, 
                'timestep': 2,  #unit in fs, could be 4 fs, but it may cause unstablity of MD
                'nsteps': 100,
                'niterations': 5000,
                'current_group_nb': None,
                'current_group_chg': None,
                'pdbx': None,
                'save_traj': False,
                'reportstate_freq': 1000,
                'savetraj_freq': 5000,
                'input_state': None,
                'if_min_heat_density': False,
                'annihilate_electrostatics': True, # if set to False, this will decouple electrostatics instead of annihilate
                'cal_adj_ene_num': 5, # To set the number of the adjacent windows of current windows to calculate the internal energy, eg. cal_adj_ene_num=5, normally will be 11 windows energy calculated, if cal_adj_ene_num=all, calculate all windows energy.
                'kbond': 10,
                'ktheta1': 10,
                'ktheta2': 10,
                'kphi1':10,
                'kphi2':10,
                'kphi3':10,
                'alchemical_co_ion': False, # should be False or int, the residue number of the co-ion (0-based), if set to False, will not add the co-ion
                },
            'restraint': {
                'restraint': False, 
                'temperature': 298.0, 
                'timestep': 2,  #unit in fs, could be 4 fs, but it may cause unstablity of MD
                'ligand_resname': 'MOL',
                'iflog_restraint_detail': False,
                'heat_nsteps' : 25000 ,
                'heat_iterations' : 50 ,
                'density_nsteps' : 25000, # density
                'npt_nsteps' : 500000, # production with plumed
                'f_npt_state' : 'npt_final_state.xml', # npt_state file
                'lambdas_json': None,
                'lambdas_group': None,
                'fake_state_xml': 'state_s0.xml', # fake state xml file for passing the simulation of state_s0
                'first_state_csv': 'state_s0.csv', # the ene loging file of free state
                'save_traj' : True,
                'reportstate_freq': 1000,
                'savetraj_freq': 5000,
                'f_plumed_input': 'plumed.dat',
                'f_plumed_output': 'Colvar',
                'plumed_freq': 100 ,# plumed record frequency
                'f_restraint': 'res_databystd.csv',
                'f_restraint2': 'res_databystd2.csv',
                'res_sele_strategy': 'lig_shape|HB_pair|HB_mainchain|Huggins',
                'fix_lig_3atoms_index': None, # '2126|2149|2147' should be used with 'res_sele_strategy' = 'Fix_lig_3atoms'
                'opt_cost_name': False, # could be 'RED_E_cost' or 'dG_forward' or False
                'if_mean': True, # if_mean is True, use the mean values of six parameters as the eq values for the harmonic restraint
                'if_init_pose': False,  # if_mean is True, use the initial values of six parameters from the input pose as the eq values for the harmonic restraint
                'preliminary_md_inital_state': None,
                'preliminary_min_and_heat': True,
                'crd': None, # assigned by the openmm-FEP-run.py -c option, but you can give the coordinate of the structure that you like for the generation of candidate groups of restraint
                'top': None, # assigned by the openmm-FEP-run.py -p option, but you can give the topology of the structure that you like for the generation of candidate groups of restraint
                }
        }


if __name__ == '__main__':
    def main():
        filename = 'input.txt'
        parser = InputParser(filename)

        # Test the md_control section
        md_control = parser.get_md_control()
        print('md_control section:', md_control)

        # Test the md_system section
        md_system = parser.get_md_system()

        # Test the pme section
        pme = parser.get_pme()


    main()

