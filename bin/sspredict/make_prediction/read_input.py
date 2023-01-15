# os 
# json
import os
import json
import pandas as pd
import numpy as np

class master_read_input:

    def __init__(self,fh):
        self.file = fh
        self.data = json.load(open(fh))
        try:
            self.model = self.data["model"]["name"]
        except:
            self.model = self.data["model"]
    def assign_model_mode(self):

        models_FCC_BCC_edge = ['FCC_Varvenne-Curtin-2016','BCC_edge_Maresca-Curtin-2019']
        model_BCC_screw_Curtin = ['BCC_screw_Maresca-Curtin-2019']
        model_BCC_screw_Suzuki = ['BCC_screw_Suzuki_RWASM-2020']
        if self.model in models_FCC_BCC_edge:
            if "pseudo-ternary" in self.data: 
                self.mode = 'edge_ternary'
            elif "compositions" in self.data:
                self.mode = 'edge_single'
            else:
                print('Please indicate correct model and calculation type in the input file.\n'
                      'Use --format to check sample inputs.')
        elif self.model in model_BCC_screw_Curtin:
            if "pseudo-ternary" in self.data: 
                self.mode = 'screw_ternary'
            elif "properties" in self.data:
                self.mode = 'screw_single'
            else:
                print('Please indicate correct model and calculation type in the input file.\n'
                      'Use --format to check sample inputs.')
        elif self.model in model_BCC_screw_Suzuki: 
            if "pseudo-ternary" in self.data:
                self.mode = "screw_suzuki_ternary"
            elif "c" in self.data['elements'][list(self.data['elements'].keys())[0]]:
                self.mode = "screw_suzuki_single"
            else:
                print('Please indicate correct model and calculation type in the input file.\n'
                      'Use --format to check sample inputs.')


class read_inputjson_edge_pseudo_ternary:
    
    def __init__(self,fh):
        if os.path.isfile(fh):
            self.file = fh
            self.data = json.load(open(fh))
            self.model = self.data["model"]
            self.name = self.data["material"]
            
        else:
            print('input file not in current directory')
            quit()

    def check_integrity_curtin_edge(self):
        def common_elements(list_1, list_2):
            a_set = set(list_1)
            b_set = set(list_2)
            if (a_set & b_set):
                return True
            else:
                return False
        def check_input_data(element_list,elements_data):
            for element in element_list:
                if not elements_data[element]:
                    print('No data found: {}'.format(element))
                
        if not self.data["pseudo-ternary"]:
            print("please supply valid pseudo-ternary")
        elif not self.data["pseudo-ternary"]["psA"]:
            print("please supply valid pseudo-ternary input: psA")
        elif not self.data["pseudo-ternary"]["psB"]:
            print("please supply valid pseudo-ternary input: psB")
        elif not self.data["pseudo-ternary"]["psC"]:
            print("please supply valid pseudo-ternary input: psC")
        if common_elements(self.data['pseudo-ternary']['psA']['grouped_elements'],self.data['pseudo-ternary']['psB']['grouped_elements']):
            print('Same element found in two groups')
        elif common_elements(self.data['pseudo-ternary']['psB']['grouped_elements'],self.data['pseudo-ternary']['psC']['grouped_elements']):
            print('Same element found in two groups')
        elif common_elements(self.data['pseudo-ternary']['psA']['grouped_elements'],self.data['pseudo-ternary']['psC']['grouped_elements']):
            print('Same element found in two groups')
        
        check_input_data(self.data['pseudo-ternary']['psA']['grouped_elements'],self.data['elements'])
        check_input_data(self.data['pseudo-ternary']['psB']['grouped_elements'],self.data['elements'])
        check_input_data(self.data['pseudo-ternary']['psC']['grouped_elements'],self.data['elements'])
    
    def grab_properties_curtin_edge(self):
        self.psA_elements = self.data['pseudo-ternary']['psA']['grouped_elements']
        self.psB_elements = self.data['pseudo-ternary']['psB']['grouped_elements']
        self.psC_elements = self.data['pseudo-ternary']['psC']['grouped_elements']
        
        self.psA_ratio = self.data['pseudo-ternary']['psA']['ratio']
        self.psB_ratio = self.data['pseudo-ternary']['psB']['ratio']
        self.psC_ratio = self.data['pseudo-ternary']['psC']['ratio']
        
        self.psA_range = self.data['pseudo-ternary']['psA']['range']
        self.psB_range = self.data['pseudo-ternary']['psB']['range']
        self.psC_range = self.data['pseudo-ternary']['psC']['range']
        
        self.increment = self.data['pseudo-ternary']['increment']
        
        self.elements_data = self.data['elements']
        
        self.elements_ABC = self.psA_elements + self.psB_elements + self.psC_elements
        
        try:
            self.temperature = np.arange(self.data['conditions']['temperature']['min'],
                                        self.data['conditions']['temperature']['max']+self.data['conditions']['temperature']['inc'],
                                        self.data['conditions']['temperature']['inc'])
        except:
            self.temperature = np.array([float(self.data['conditions']['temperature'])])
            
        self.strain_r = self.data['conditions']['strain_r']
        self.exp_conditions = [self.temperature,self.strain_r]
        self.structure = (self.data["structure"].lower())
        self.model = self.data['model']['name']
        if self.model in ['FCC_Varvenne-Curtin-2016','BCC_edge_Maresca-Curtin-2019']:
            if 'f1' in self.data['model']['name']:
                self.f1 = self.data['model']['name']['f1']
            else:
                if self.structure == 'fcc':
                    self.f1 = 0.35
                else: 
                    self.f1 = 0.7844
            if 'f2' in self.data['model']['name']:
                self.f2 = self.data['model']['name']['f2']
            else:
                if self.structure == 'fcc':
                    self.f2 = 5.7
                else: 
                    self.f2 = 7.2993
                
            if 'alpha' in self.data['model']['name']:
                self.alpha = self.data['model']['name']['alpha']
            else:
                self.alpha = 0.123
        self.dislocation_properties = [self.alpha,self.f1,self.f2]
        for element_i in self.elements_ABC:
            # compute E, nu, G for elements if not supplied
            # two of the E/nu/G must be supplied to calculate the missing one
            if not 'nu' in self.elements_data[element_i]: 
                self.elements_data[element_i]['nu'] = round(self.elements_data[element_i]['E']/2/self.elements_data[element_i]['G'] - 1,3)
            if not 'G' in self.elements_data[element_i]:
                self.elements_data[element_i]['G'] = round(self.elements_data[element_i]['E']/2/(self.elements_data[element_i]['nu'] + 1),1)
            if not 'E' in self.elements_data[element_i]:
                self.elements_data[element_i]['E'] = round(self.elements_data[element_i]['G']*2*(self.elements_data[element_i]['nu'] + 1),1)
            
            # compute a/b/Vn for elements based on lattice structure
            # one of the a/b/Vn must be supplied to compute the other two
            if self.structure == 'fcc':
                # fcc: Vn = a^3/4, b = a/sqrt(2)
                if 'Vn' in self.elements_data[element_i]:
                    self.elements_data[element_i]['a'] = (self.elements_data[element_i]['Vn']*4)**(1/3)
                    self.elements_data[element_i]['b'] = self.elements_data[element_i]['a']/np.sqrt(2)
                elif 'b' in self.elements_data[element_i]:
                    self.elements_data[element_i]['a'] = self.elements_data[element_i]['b']*np.sqrt(2)
                    self.elements_data[element_i]['Vn'] = self.elements_data[element_i]['a']**3/4
                elif 'a' in self.elements_data[element_i]:
                    self.elements_data[element_i]['b'] = self.elements_data[element_i]['a']/np.sqrt(2)
                    self.elements_data[element_i]['Vn'] = self.elements_data[element_i]['a']**3/4
            else: 
                # bcc: Vn = a^3/2, b = a*sqrt(3)/2
                if 'Vn' in self.elements_data[element_i]:
                    self.elements_data[element_i]['a'] = (self.elements_data[element_i]['Vn']*2)**(1/3)
                    self.elements_data[element_i]['b'] = self.elements_data[element_i]['a']*np.sqrt(3)/2
                elif 'b' in self.elements_data[element_i]:
                    self.elements_data[element_i]['a'] = self.elements_data[element_i]['b']*2/np.sqrt(3)
                    self.elements_data[element_i]['Vn'] = self.elements_data[element_i]['a']**3/2
                elif 'a' in self.elements_data[element_i]:
                    self.elements_data[element_i]['b'] = self.elements_data[element_i]['a']*np.sqrt(3)/2
                    self.elements_data[element_i]['Vn'] = self.elements_data[element_i]['a']**3/2
        if "uncertainty_level" in self.data:
            #'a': uncertainty of lattice constant, if on, default is 1% 
            #'elastic constant': uncertainty of lattice constants, if on, default is 5% for each element
            
            if self.data["uncertainty_level"]["on/off"].lower() == "on":
                if "a" in self.data["uncertainty_level"]:
                    self.uncertainty_a = self.data["uncertainty_level"]['a']
                else: 
                    self.uncertainty_a = 0.01
                if "elastic_constants" in self.data["uncertainty_level"]:
                    self.uncertainty_EGv = self.data["uncertainty_level"]['elastic_constants']
                else: 
                    self.uncertainty_EGv = 0.05
            else: 
                self.uncertainty_a = 0
                self.uncertainty_EGv = 0
        self.uncertainty_levels = [self.uncertainty_a,self.uncertainty_EGv]
        try: 
            self.savefilename = self.data['savefile']
        except:
            self.savefilename = self.data["material"] + '_out'
        


class read_inputjson_edge_single_calculation:
# group the elements from reading the input file
    def __init__(self,fh):
        if os.path.isfile(fh):
            self.file = fh
            self.data = json.load(open(fh))
            self.model = self.data["model"]
            self.name = self.data["material"]
        else:
            print('input file not in current directory')
            quit()

    def check_integrity_curtin_edge(self):
        def common_elements(list_1, list_2):
            a_set = set(list_1)
            b_set = set(list_2)
            if (a_set & b_set):
                return True
            else:
                return False
        def check_input_data(element_list,elements_data):
            for element in element_list:
                if not elements_data[element]:
                    print('No data found: {}'.format(element))
                
    def grab_properties_curtin_edge(self):

        
        self.elements_data = self.data['elements']
        
        self.elements_order = self.data['compositions']['element_order']
            
        self.concentrations = np.array(self.data['compositions']['concentrations'])
        
        try:
            self.temperature = np.arange(self.data['conditions']['temperature']['min'],
                                        self.data['conditions']['temperature']['max']+self.data['conditions']['temperature']['inc'],
                                        self.data['conditions']['temperature']['inc'])
        except:
            self.temperature = np.array([float(self.data['conditions']['temperature'])])
            
        element_composition = {}
        for i in range(len(self.elements_order)):
            element_composition[self.elements_order[i]] = self.concentrations.T[i]
        self.element_composition = pd.DataFrame(data=element_composition)
        self.strain_r = self.data['conditions']['strain_r']
        self.exp_conditions = [self.temperature,self.strain_r]
        self.structure = (self.data["structure"].lower())
        self.model = self.data['model']['name']
        if self.model in ['FCC_Varvenne-Curtin-2016','BCC_edge_Maresca-Curtin-2019']:
            if 'f1' in self.data['model']:
                self.f1 = self.data['model']['f1']
            else:
                if self.structure == 'fcc':
                    self.f1 = 0.35
                else: 
                    self.f1 = 0.7844
            if 'f2' in self.data['model']:
                self.f2 = self.data['model']['f2']
            else:
                if self.structure == 'fcc':
                    self.f2 = 5.7
                else: 
                    self.f2 = 7.2993
                
            if 'alpha' in self.data['model']:
                self.alpha = self.data['model']['alpha']
            else:
                self.alpha = 0.123
        self.adjustable_paras= [self.alpha,self.f1,self.f2]
        for element_i in self.elements_order:
            # compute E, nu, G for elements if not supplied
            # two of the E/nu/G must be supplied to calculate the missing one
            if not 'nu' in self.elements_data[element_i]: 
                self.elements_data[element_i]['nu'] = round(self.elements_data[element_i]['E']/2/self.elements_data[element_i]['G'] - 1,3)
            if not 'G' in self.elements_data[element_i]:
                self.elements_data[element_i]['G'] = round(self.elements_data[element_i]['E']/2/(self.elements_data[element_i]['nu'] + 1),1)
            if not 'E' in self.elements_data[element_i]:
                self.elements_data[element_i]['E'] = round(self.elements_data[element_i]['G']*2*(self.elements_data[element_i]['nu'] + 1),1)
            
            # compute a/b/Vn for elements based on lattice structure
            # one of the a/b/Vn must be supplied to compute the other two
            if self.structure == 'fcc':
                # fcc: Vn = a^3/4, b = a/sqrt(2)
                if 'Vn' in self.elements_data[element_i]:
                    self.elements_data[element_i]['a'] = (self.elements_data[element_i]['Vn']*4)**(1/3)
                    self.elements_data[element_i]['b'] = self.elements_data[element_i]['a']/np.sqrt(2)
                elif 'b' in self.elements_data[element_i]:
                    self.elements_data[element_i]['a'] = self.elements_data[element_i]['b']*np.sqrt(2)
                    self.elements_data[element_i]['Vn'] = self.elements_data[element_i]['a']**3/4
                elif 'a' in self.elements_data[element_i]:
                    self.elements_data[element_i]['b'] = self.elements_data[element_i]['a']/np.sqrt(2)
                    self.elements_data[element_i]['Vn'] = self.elements_data[element_i]['a']**3/4
            else: 
                # bcc: Vn = a^3/2, b = a*sqrt(3)/2
                if 'Vn' in self.elements_data[element_i]:
                    self.elements_data[element_i]['a'] = (self.elements_data[element_i]['Vn']*2)**(1/3)
                    self.elements_data[element_i]['b'] = self.elements_data[element_i]['a']*np.sqrt(3)/2
                elif 'b' in self.elements_data[element_i]:
                    self.elements_data[element_i]['a'] = self.elements_data[element_i]['b']*2/np.sqrt(3)
                    self.elements_data[element_i]['Vn'] = self.elements_data[element_i]['a']**3/2
                elif 'a' in self.elements_data[element_i]:
                    self.elements_data[element_i]['b'] = self.elements_data[element_i]['a']*np.sqrt(3)/2
                    self.elements_data[element_i]['Vn'] = self.elements_data[element_i]['a']**3/2
                    
        if "uncertainty_level" in self.data:
            #'a': uncertainty of lattice constant, if on, default is 1% 
            #'elastic constant': uncertainty of lattice constants, if on, default is 5% for each element
            
            if self.data["uncertainty_level"]["on/off"].lower() == "on":
                if "a" in self.data["uncertainty_level"]:
                    self.uncertainty_a = self.data["uncertainty_level"]['a']
                else: 
                    self.uncertainty_a = 0.01
                if "elastic_constants" in self.data["uncertainty_level"]:
                    self.uncertainty_EGv = self.data["uncertainty_level"]['elastic_constants']
                else: 
                    self.uncertainty_EGv = 0.05
            else: 
                self.uncertainty_a = 0
                self.uncertainty_EGv = 0
        self.uncertainty_levels = [self.uncertainty_a,self.uncertainty_EGv]
    


class read_inputjson_BCC_screw_pseudo_ternary: 
    def __init__(self,fh):
        if os.path.isfile(fh):
            self.file = fh
            self.data = json.load(open(fh))
        else:
            print('input file not in current directory')
            quit()
        
        self.model = self.data["model"]
        self.name = self.data["material"]
        # properties
        self.elements_data = self.data['elements']
        # pseudo-ternary
        self.psA_elements = self.data['pseudo-ternary']['psA']['grouped_elements']
        self.psB_elements = self.data['pseudo-ternary']['psB']['grouped_elements']
        self.psC_elements = self.data['pseudo-ternary']['psC']['grouped_elements']
        
        self.psA_ratio = self.data['pseudo-ternary']['psA']['ratio']
        self.psB_ratio = self.data['pseudo-ternary']['psB']['ratio']
        self.psC_ratio = self.data['pseudo-ternary']['psC']['ratio']
        
        self.psA_range = self.data['pseudo-ternary']['psA']['range']
        self.psB_range = self.data['pseudo-ternary']['psB']['range']
        self.psC_range = self.data['pseudo-ternary']['psC']['range']
        self.increment = self.data['pseudo-ternary']['increment']
        
        # adjustable paras
        self.adjustable_paras = self.data['adjustables']
        # exp conditions
        self.conditions = self.data['conditions']

        # output file
        try: 
            self.savefilename = self.data['savefile']
        except:
            self.savefilename = self.data["material"] + '_out'


class read_inputjson_BCC_screw_single_calculation: 
    def __init__(self,fh):
        if os.path.isfile(fh):
            self.file = fh
            self.data = json.load(open(fh))
        else:
            print('input file not in current directory')
            quit()
        
        self.model = self.data["model"]
        self.name = self.data["material"]
        # properties
        self.properties = self.data['properties']
        # adjustable paras
        self.adjustable_paras = self.data['adjustables']
        # exp conditions
        self.conditions = self.data['conditions']
        # output file
        try: 
            self.savefilename = self.data['savefile']
        except:
            self.savefilename = self.data["material"] + '_out'



class read_json_Suzuki_model_RWASM_ternary:
    
    def __init__(self,fh):
        if os.path.isfile(fh):
            self.file = fh
            self.data = json.load(open(fh))
        else:
            print('input file not in current directory')
            quit()
        
        self.model = self.data["model"]
        self.name = self.data["material"]
        # properties
        self.elements_data = self.data['elements']
        # adjustable paras
        self.adjustable_paras = self.data['adjustables']
        # exp conditions
        self.conditions = self.data['conditions']
        # output file
        try: 
            self.savefilename = self.data['savefile']
        except:
            self.savefilename = self.data["material"] + '_out'
        # ternary 
        self.psA_elements = self.data['pseudo-ternary']['psA']['grouped_elements']
        self.psB_elements = self.data['pseudo-ternary']['psB']['grouped_elements']
        self.psC_elements = self.data['pseudo-ternary']['psC']['grouped_elements']
        
        self.psA_ratio = self.data['pseudo-ternary']['psA']['ratio']
        self.psB_ratio = self.data['pseudo-ternary']['psB']['ratio']
        self.psC_ratio = self.data['pseudo-ternary']['psC']['ratio']
        
        self.psA_range = self.data['pseudo-ternary']['psA']['range']
        self.psB_range = self.data['pseudo-ternary']['psB']['range']
        self.psC_range = self.data['pseudo-ternary']['psC']['range']
        
        self.increment = self.data['pseudo-ternary']['increment']
        # check if liquidus tempratures are supplied
        try:
            self.T_l_tag = self.data['TC_data']['T_liquidus']
            if self.T_l_tag == 'tc_liquidus_wrapper':
                print('Will calculate liquidus temperatures from TC-Python.')
                print('Make sure to run with desired TC database. use TC_database tag under TC_data. ')
                try: 
                    self.tcdatabase = self.data['TC_data']['TC_database']
                    print('current database: {}'.format(self.data['TC_data']['TC_database']))
                except:
                    self.tcdatabase = 'TCHEA4'
                    print('current database: TCHEA4 (default)')

            else:
                try:
                    T_l_df = pd.read_csv(self.T_l_tag)
                    self.T_l = T_l_df['T_liquidus']
                    print('Getting liquidus temperatures from file {}'.format(self.T_l_tag))
                except: 
                    self.T_l = None
                    print('Cannot read the file {}, be sure it is .csv format containing the \'T_liquidus\' '.format(self.T_l_tag))
        except:
            self.T_l = None
            self.T_l_tag = None


        
class read_json_Suzuki_model_RWASM_T:
    
    def __init__(self,fh):
        if os.path.isfile(fh):
            self.file = fh
            self.data = json.load(open(fh))
        else:
            print('input file not in current directory')
            quit()
        
        self.model = self.data["model"]
        self.name = self.data["material"]
        # properties
        self.elements_data = self.data['elements']
        # compositions
        self.elements_order = self.data['compositions']['element_order']
        self.concentrations = self.data['compositions']['concentrations']
        #
        conditions = self.data['conditions']
        self.strain_r = conditions['strain_r']
        self.temperature = np.arange(conditions['temperature']['min'],
                               conditions['temperature']['max']+conditions['temperature']['inc'],
                               conditions['temperature']['inc'])
        # exp conditions
        self.conditions = [self.temperature,self.strain_r]
        
        element_composition = {}
        for i in range(len(self.elements_order)):
            c_x = np.array(self.concentrations).transpose()[i]
            c_x_T = np.array(c_x).transpose().flatten()
            element_composition[self.elements_order[i]] = c_x_T
        self.element_composition = pd.DataFrame(data=element_composition)
        self.T_l = np.ones(len(self.element_composition))*1e8
        try:
            T_l = np.genfromtxt(np.array(self.data['jog_drag_inputs']['T_l']))
        except:
            T_l = None
            print('No liquidus temperature supplied. Not perform high-T jog-dragging calculations.')
        if T_l is not None:
            for i in range(len(self.element_composition)):
                try:
                    self.T_l[i]=T_l[i]
                except:
                    print('Not perform high-T jog-dragging calculation of {}={} because no liquidus temperature supplied.'.format(self.element_composition.columns.values,self.element_composition.iloc[i].values))
        self.T_l = np.nan_to_num(self.T_l,nan=1e8)
        # adjustable paras
        self.adjustable_paras = self.data['adjustables']
        
        # output file
        try: 
            self.savefilename = self.data['savefile']
        except:
            self.savefilename = self.data["material"] + '_out'
        


class make_adjustables:


    def __init__(self,models,
                        f1=None,f2=None,alpha=None, # edge models
                        kink_width=None,Delta_V_p_para=None,Delta_E_p_para=None,# MC screw model
                        tau_i_exponent=None,dislocation_density=None,trial_kappa=None,trial_tau_k=None # Suzuki screw model
                        ):

        # adjustables for different models
        self.adjustables = {}
        self.ssmodels = models # list of models 
        # if using "FCC_Varvenne-Curtin-2016"
        if "FCC_Varvenne-Curtin-2016" in self.ssmodels:
            print('Setting adjustable parameters of the model {}.'.format("FCC_Varvenne-Curtin-2016"))
            self.adjustables_VC = {
                "f1":   0.35,  # dimensionless pressure field parameter for athermal yield stress 
                "f2":   5.7,   # dimensionless pressure field parameter for energy barrier
                "alpha": 0.123 # dislocation line tension parameter
                }
            if f1 is not None:
                self.adjustables_VC['f1'] = float(f1)
            if f2 is not None:
                self.adjustables_VC['f2'] = float(f2)
            if alpha is not None:
                self.adjustables_VC['alpha'] = float(alpha)
            print(self.adjustables_VC)
            self.adjustables["FCC_Varvenne-Curtin-2016"] = self.adjustables_VC
        # if using "BCC_edge_Maresca-Curtin-2019"
        if "BCC_edge_Maresca-Curtin-2019" in self.ssmodels:
            print('Setting adjustable parameters of the model {}.'.format("BCC_edge_Maresca-Curtin-2019"))
            self.adjustables_MC_edge = {
                "f1":   0.7844,# dimensionless pressure field parameter for athermal yield stress 
                "f2":   7.2993,# dimensionless pressure field parameter for energy barrier
                "alpha": 0.123 # dislocation line tension parameter
                }
            if f1 is not None:
                self.adjustables_MC_edge['f1'] = float(f1)
            if f2 is not None:
                self.adjustables_MC_edge['f2'] = float(f2)
            if alpha is not None:
                self.adjustables_MC_edge['alpha'] = float(alpha)
            print(self.adjustables_MC_edge)
            self.adjustables["BCC_edge_Maresca-Curtin-2019"] = self.adjustables_MC_edge
        # if using "BCC_screw_Maresca-Curtin-2019"
        if "BCC_screw_Maresca-Curtin-2019" in self.ssmodels:
            print('Setting adjustable parameters of the model {}.'.format("BCC_screw_Maresca-Curtin-2019"))
            self.adjustables_MC_screw = {
                "kink_width":10,      # kink width para, usually 10b to 20b, b=burgers vector
                "Delta_V_p_para":1, # Peierls barrier para
                "Delta_E_p_para":1  # solute-dislocation interaction energy para
                }
            if kink_width is not None:
                self.adjustables_MC_screw['kink_width'] = float(kink_width)
            if Delta_V_p_para is not None:
                self.adjustables_MC_screw['Delta_V_p_para'] = float(Delta_V_p_para)
            if Delta_E_p_para is not None:
                self.adjustables_MC_screw['Delta_E_p_para'] = float(Delta_E_p_para)
            print(self.adjustables_MC_screw)
            self.adjustables["BCC_screw_Maresca-Curtin-2019"] = self.adjustables_MC_screw
        # if using "BCC_screw_Suzuki_RWASM-2020"
        if "BCC_screw_Suzuki_RWASM-2020" in self.ssmodels:
            print('Setting adjustable parameters of the model {}.'.format("BCC_screw_Suzuki_RWASM-2020"))
            self.adjustables_SRWASM_screw = {
                    "kink_width":10,        # kink width para, usually 10b to 20b, b=burgers vector
                    "tau_i_exponent":1,     # exponent for phenomenological summation of yield stresses of components, usually 0.9-1
                    "dislocation_density":4e13, # dislocation density 
                    "trial_kappa":{         # trial kappa range for optimization 
                        "min":1,            
                        "max":4,
                        "inc":0.5         
                    },
                    "trial_tau_k":5         # trial tau_k for root finding, unit: MPa
                }
            if kink_width is not None:
                self.adjustables_SRWASM_screw['kink_width'] = float(kink_width)
            if tau_i_exponent is not None:
                self.adjustables_SRWASM_screw['tau_i_exponent'] = float(tau_i_exponent)
            if trial_tau_k is not None:
                self.adjustables_SRWASM_screw['trial_tau_k'] = float(trial_tau_k)
            if trial_kappa is not None:
                try:
                    self.adjustables_SRWASM_screw['trial_kappa']['min'] = float(trial_kappa['min'])
                    self.adjustables_SRWASM_screw['trial_kappa']['max'] = float(trial_kappa['max'])
                    self.adjustables_SRWASM_screw['trial_kappa']['inc'] = float(trial_kappa['inc'])
                except:
                    print('trial_kappa need to be a dictionary containing min, max, and inc. e.g. trial_kappa={"min":1,"max":4,"inc":0.5}')
            print(self.adjustables_SRWASM_screw)
            self.adjustables["BCC_screw_Suzuki_RWASM-2020"] = self.adjustables_SRWASM_screw
        

class get_elements_data:

    def __init__(self,elements,models,fh=None,from_dict=None,alias=None,averaging_scheme='default'):
        self.ssmodels = models 
        # read from file or from dictionary 
        if fh is not None: 
            print('reading datafile {}.'.format(fh))
            try:
                js = json.load(open(fh))
                self.data = js['elements']
            except:
                print('Please provide a valid JSON.')
        elif from_dict is not None:
            try:
                self.data = from_dict['elements']
            except:
                print('Please provide a valid dictionary.')
        else:
            print('Please provide a valid JSON or dictionary.')
        #print(self.data)
        # setting up data for each selected model
        self.data_of_ssmodels = {}
        self.averaging_scheme=averaging_scheme 
        print('grab elemental data for prediction.')
        
        if "FCC_Varvenne-Curtin-2016" in self.ssmodels:
            print('Setting elemental data of the model {}.'.format("FCC_Varvenne-Curtin-2016"))
            # crystal data: Vn/a/b, structure
            # elastic constants 
            elements_data = {}
            try:
                for element_i in elements: 
                    elements_data[element_i] = {}
                    if alias is not None:
                        try:
                            element_i_alias = alias[element_i]
                        except: 
                            element_i_alias = element_i
                    else:
                        element_i_alias = element_i
                    # first get elastic constants, 
                    # two of the E/nu/G must be supplied to compute the rest one
                    #### or supply the stiffness matrix + averaging scheme ----2022.10.21
                    if self.averaging_scheme == 'default': 
                        elements_data['averaging_scheme'] = self.averaging_scheme
                        try:
                            elements_data[element_i]['nu'] = round(self.data[element_i_alias]['nu'],3)
                        except:
                            elements_data[element_i]['nu'] = round(self.data[element_i_alias]['E']/2/self.data[element_i_alias]['G'] - 1,3)
                        try:
                            elements_data[element_i]['G'] = round(self.data[element_i_alias]['G'],1)
                        except:
                            elements_data[element_i]['G'] = round(self.data[element_i_alias]['E']/2/(self.data[element_i_alias]['nu'] + 1),1)
                        try:
                            elements_data[element_i]['E'] = round(self.data[element_i_alias]['E'],1)
                        except:
                            elements_data[element_i]['E'] = round(self.data[element_i_alias]['G']*2*(self.data[element_i_alias]['nu'] + 1),1)
                    else: 
                        # get the stiffness matrix 6x6
                        elements_data[element_i]['Cij'] = np.round(np.array(self.data[element_i_alias]['Cij']),2)
                        elements_data['averaging_scheme'] = self.averaging_scheme

                    if self.data[element_i_alias]['structure'] == 'fcc':
                        # fcc: Vn = a^3/4, b = a/sqrt(2)
                        if 'Vn' in self.data[element_i_alias]:
                            elements_data[element_i]['Vn'] = self.data[element_i_alias]['Vn']
                        if 'b' in self.data[element_i_alias]:
                            elements_data[element_i]['b'] = self.data[element_i_alias]['b']
                        if 'a' in self.data[element_i_alias]:
                            elements_data[element_i]['a'] = self.data[element_i_alias]['a']

                        if "Vn" not in elements_data[element_i]:
                            try: 
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                            except:
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*np.sqrt(2)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                        if "a" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                            except: 
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*np.sqrt(2)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                        if "b" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)
                            except:
                                elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                                elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)
                    else: 
                        # bcc: Vn = a^3/2, b = a*sqrt(3)/2

                        if 'Vn' in self.data[element_i_alias]:
                            elements_data[element_i]['Vn'] = self.data[element_i_alias]['Vn']
                        if 'b' in self.data[element_i_alias]:
                            elements_data[element_i]['b'] = self.data[element_i_alias]['b']
                        if 'a' in self.data[element_i_alias]:
                            elements_data[element_i]['a'] = self.data[element_i_alias]['a']

                        if "Vn" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                            except:
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*2/np.sqrt(3)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                        if "a" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['a'] = np.cbrt(elements_data[element_i]['Vn']*2)
                            except: 
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*2/np.sqrt(3)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                        if "b" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                            except:
                                elements_data[element_i]['a'] = np.cbrt(elements_data[element_i]['Vn']*2)
                                elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                        # because using FCC model, lattice constant and b should be in the FCC medium: 
                        # the atomic volume Vn is considered unchanged for a BCC element in FCC structure 
                        elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                        elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)

            except:
                print('Failed to fetch data for {}'.format(elements))
                print('Make sure data file contains: ')
                print('1. "structure", \n'
                      '2. "Vn" or "a" or "b",  \n '
                      '3. "E" "G" "nu"  or "Cij" \n ')
            self.data_of_ssmodels["FCC_Varvenne-Curtin-2016"] = elements_data
        # if using "BCC_edge_Maresca-Curtin-2019"
        if "BCC_edge_Maresca-Curtin-2019" in self.ssmodels:
            print('Setting elemental data of the model {}.'.format("BCC_edge_Maresca-Curtin-2019"))
            # crystal data: Vn/a/b, structure
            # elastic constants
            elements_data = {}
            try:
                for element_i in elements: 
                    if alias is not None:
                        try:
                            element_i_alias = alias[element_i]
                        except: 
                            element_i_alias = element_i
                    else:
                        element_i_alias = element_i
                    elements_data[element_i] = {}
                    # first get elastic constants, 
                    # two of the E/nu/G must be supplied to compute the rest one
                    #### or supply the stiffness matrix + averaging scheme ----2022.10.21
                    if self.averaging_scheme == 'default': 
                        elements_data['averaging_scheme'] = self.averaging_scheme
                        try:
                            elements_data[element_i]['nu'] = round(self.data[element_i_alias]['nu'],3)
                        except:
                            elements_data[element_i]['nu'] = round(self.data[element_i_alias]['E']/2/self.data[element_i_alias]['G'] - 1,3)
                        try:
                            elements_data[element_i]['G'] = round(self.data[element_i_alias]['G'],1)
                        except:
                            elements_data[element_i]['G'] = round(self.data[element_i_alias]['E']/2/(self.data[element_i_alias]['nu'] + 1),1)
                        try:
                            elements_data[element_i]['E'] = round(self.data[element_i_alias]['E'],1)
                        except:
                            elements_data[element_i]['E'] = round(self.data[element_i_alias]['G']*2*(self.data[element_i_alias]['nu'] + 1),1)
                    else: 
                        # get the stiffness matrix 6x6
                        elements_data[element_i]['Cij'] = np.round(np.array(self.data[element_i_alias]['Cij']),2)
                        elements_data['averaging_scheme'] = self.averaging_scheme

                    if self.data[element_i_alias]['structure'] == 'fcc':
                        # fcc: Vn = a^3/4, b = a/sqrt(2)
                        if 'Vn' in self.data[element_i_alias]:
                            elements_data[element_i]['Vn'] = self.data[element_i_alias]['Vn']
                        if 'b' in self.data[element_i_alias]:
                            elements_data[element_i]['b'] = self.data[element_i_alias]['b']
                        if 'a' in self.data[element_i_alias]:
                            elements_data[element_i]['a'] = self.data[element_i_alias]['a']

                        if "Vn" not in elements_data[element_i]:
                            try: 
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                            except:
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*np.sqrt(2)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                        if "a" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                            except: 
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*np.sqrt(2)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                        if "b" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)
                            except:
                                elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                                elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)
                        # because using BCC model, lattice constant and b should be in the BCC medium: 
                        # the atomic volume Vn is considered unchanged for a FCC element in BCC structure 
                        elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*2)**(1/3)
                        elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                    else: 
                        # bcc: Vn = a^3/2, b = a*sqrt(3)/2

                        if 'Vn' in self.data[element_i_alias]:
                            elements_data[element_i]['Vn'] = self.data[element_i_alias]['Vn']
                        if 'b' in self.data[element_i_alias]:
                            elements_data[element_i]['b'] = self.data[element_i_alias]['b']
                        if 'a' in self.data[element_i_alias]:
                            elements_data[element_i]['a'] = self.data[element_i_alias]['a']

                        if "Vn" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                            except:
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*2/np.sqrt(3)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                        if "a" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['a'] = np.cbrt(elements_data[element_i]['Vn']*2)
                            except: 
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*2/np.sqrt(3)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                        if "b" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                            except:
                                elements_data[element_i]['a'] = np.cbrt(elements_data[element_i]['Vn']*2)
                                elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2

            except:
                print('Failed to fetch data for {}'.format(elements))
                print('Make sure data file contains: ')
                print('1. "structure", \n'
                      '2. "Vn" or "a" or "b",  \n '
                      '3. "E" "G" "nu" or "Cij" \n ')
            self.data_of_ssmodels["BCC_edge_Maresca-Curtin-2019"] = elements_data
        # if using "BCC_screw_Maresca-Curtin-2019"
        if "BCC_screw_Maresca-Curtin-2019" in self.ssmodels:
            print('Setting elemental data of the model {}.'.format("BCC_screw_Maresca-Curtin-2019"))
            # 1. crystal data: Vn/a/b, structure
            # 2. kink formation energy
            # 3. vacancy and self interstitial formation energy
            # 4. solute-dislocation interaction and Peierls barrier
            elements_data = {}
            try:
                for element_i in elements: 
                    elements_data[element_i] = {}
                    if alias is not None:
                        try:
                            element_i_alias = alias[element_i]
                        except: 
                            element_i_alias = element_i
                    else: 
                        element_i_alias = element_i
                    if self.data[element_i_alias]['structure'] == 'fcc':
                        # fcc: Vn = a^3/4, b = a/sqrt(2)
                        if 'Vn' in self.data[element_i_alias]:
                            elements_data[element_i]['Vn'] = self.data[element_i_alias]['Vn']
                        if 'b' in self.data[element_i_alias]:
                            elements_data[element_i]['b'] = self.data[element_i_alias]['b']
                        if 'a' in self.data[element_i_alias]:
                            elements_data[element_i]['a'] = self.data[element_i_alias]['a']

                        if "Vn" not in elements_data[element_i]:
                            try: 
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                            except:
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*np.sqrt(2)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                        if "a" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                            except: 
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*np.sqrt(2)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                        if "b" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)
                            except:
                                elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                                elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)
                        # because using BCC model, lattice constant and b should be in the BCC medium: 
                        # the atomic volume Vn is considered unchanged for a FCC element in BCC structure 
                        elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*2)**(1/3)
                        elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                    else: 
                        # bcc: Vn = a^3/2, b = a*sqrt(3)/2

                        if 'Vn' in self.data[element_i_alias]:
                            elements_data[element_i]['Vn'] = self.data[element_i_alias]['Vn']
                        if 'b' in self.data[element_i_alias]:
                            elements_data[element_i]['b'] = self.data[element_i_alias]['b']
                        if 'a' in self.data[element_i_alias]:
                            elements_data[element_i]['a'] = self.data[element_i_alias]['a']

                        if "Vn" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                            except:
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*2/np.sqrt(3)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                        if "a" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['a'] = np.cbrt(elements_data[element_i]['Vn']*2)
                            except: 
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*2/np.sqrt(3)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                        if "b" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                            except:
                                elements_data[element_i]['a'] = np.cbrt(elements_data[element_i]['Vn']*2)
                                elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                    # vacancy and self interstitial formation energy
                    elements_data[element_i]['E_f_v'] = self.data[element_i_alias]['E_f_v']
                    elements_data[element_i]['E_f_si'] = self.data[element_i_alias]['E_f_si']
                    # Delta_E_p: characteristic energy fluctuation due to solute-dislocation interaction 
                    elements_data[element_i]['Delta_E_p'] = self.data[element_i_alias]['Delta_E_p']
                    # Peierls barrier 
                    elements_data[element_i]['Delta_V_p'] = self.data[element_i_alias]['Delta_V_p']
                    elements_data[element_i]['E_k'] = self.data[element_i_alias]['E_k']

            except:
                print('Failed to fetch data for {}'.format(elements))
                print('Make sure data file contains: ')
                print('1. "structure", \n'
                      '2. "Vn" or "a" or "b",  \n '
                      '3. "E_f_v" "E_f_si" "Delta_E_p" "E_k"  \n ')
            self.data_of_ssmodels["BCC_screw_Maresca-Curtin-2019"] = elements_data
        # if using "BCC_screw_Suzuki_RWASM-2020"
        if "BCC_screw_Suzuki_RWASM-2020" in self.ssmodels:
            print('Setting elemental data of the model {}.'.format("BCC_screw_Suzuki_RWASM-2020"))
            # 1. crystal data: Vn/a/b, structure
            # 2. vacancy and self interstitial formation energy
            # 3. solute-dislocation interaction energy 
            # 4. elastic constants G and nu

            elements_data = {}
            elements_data['averaging_scheme'] = self.averaging_scheme
            try:
                for element_i in elements: 
                    if alias is not None:
                        try:
                            element_i_alias = alias[element_i]
                        except: 
                            element_i_alias = element_i
                    else:
                        element_i_alias = element_i
                    elements_data[element_i] = {}
                    # first get elastic constants, 
                    # two of the E/nu/G must be supplied to compute the rest one
                    #### or supply the stiffness matrix + averaging scheme ----2022.10.21
                    if self.averaging_scheme == 'default': 
                        #elements_data['averaging_scheme'] = self.averaging_scheme
                        try:
                            elements_data[element_i]['nu'] = round(self.data[element_i_alias]['nu'],3)
                        except:
                            elements_data[element_i]['nu'] = round(self.data[element_i_alias]['E']/2/self.data[element_i_alias]['G'] - 1,3)
                        try:
                            elements_data[element_i]['G'] = round(self.data[element_i_alias]['G'],1)
                        except:
                            elements_data[element_i]['G'] = round(self.data[element_i_alias]['E']/2/(self.data[element_i_alias]['nu'] + 1),1)
                        try:
                            elements_data[element_i]['E'] = round(self.data[element_i_alias]['E'],1)
                        except:
                            elements_data[element_i]['E'] = round(self.data[element_i_alias]['G']*2*(self.data[element_i_alias]['nu'] + 1),1)
                    else: 
                        # get the stiffness matrix 6x6
                        elements_data[element_i]['Cij'] = np.round(np.array(self.data[element_i_alias]['Cij']),2)
                        #elements_data['averaging_scheme'] = self.averaging_scheme

                    if self.data[element_i_alias]['structure'] == 'fcc':
                        # fcc: Vn = a^3/4, b = a/sqrt(2)
                        if 'Vn' in self.data[element_i_alias]:
                            elements_data[element_i]['Vn'] = self.data[element_i_alias]['Vn']
                        if 'b' in self.data[element_i_alias]:
                            elements_data[element_i]['b'] = self.data[element_i_alias]['b']
                        if 'a' in self.data[element_i_alias]:
                            elements_data[element_i]['a'] = self.data[element_i_alias]['a']

                        if "Vn" not in elements_data[element_i]:
                            try: 
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                            except:
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*np.sqrt(2)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                        if "a" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                            except: 
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*np.sqrt(2)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/4
                        if "b" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)
                            except:
                                elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*4)**(1/3)
                                elements_data[element_i]['b'] = elements_data[element_i]['a']/np.sqrt(2)
                        # because using BCC model, lattice constant and b should be in the BCC medium: 
                        # the atomic volume Vn is considered unchanged for a FCC element in BCC structure 
                        elements_data[element_i]['a'] = (elements_data[element_i]['Vn']*2)**(1/3)
                        elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                    else: 
                        # bcc: Vn = a^3/2, b = a*sqrt(3)/2

                        if 'Vn' in self.data[element_i_alias]:
                            elements_data[element_i]['Vn'] = self.data[element_i_alias]['Vn']
                        if 'b' in self.data[element_i_alias]:
                            elements_data[element_i]['b'] = self.data[element_i_alias]['b']
                        if 'a' in self.data[element_i_alias]:
                            elements_data[element_i]['a'] = self.data[element_i_alias]['a']

                        if "Vn" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                            except:
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*2/np.sqrt(3)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                        if "a" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['a'] = np.cbrt(elements_data[element_i]['Vn']*2)
                            except: 
                                elements_data[element_i]['a'] = elements_data[element_i]['b']*2/np.sqrt(3)
                                elements_data[element_i]['Vn'] = elements_data[element_i]['a']**3/2
                        if "b" not in elements_data[element_i]:
                            try:
                                elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2
                            except:
                                elements_data[element_i]['a'] = np.cbrt(elements_data[element_i]['Vn']*2)
                                elements_data[element_i]['b'] = elements_data[element_i]['a']*np.sqrt(3)/2

                    # E_w_i: solute-dislocation interaction energy 
                    elements_data[element_i]['E_w_i'] = self.data[element_i_alias]['E_w_i']
                    try:
                        # vacancy and self interstitial formation energy
                        elements_data[element_i]['E_f_v'] = self.data[element_i_alias]['E_f_v']
                        elements_data[element_i]['E_f_si'] = self.data[element_i_alias]['E_f_si']
                    except:
                        # no formation energy of defects provided, use the Frank-Read bowing for tau_j
                        elements_data[element_i]['E_f_v'] = '-'
                        elements_data[element_i]['E_f_si'] = '-'
                        print(' "E_f_v" "E_f_si" not provided, use Frank-Read bowing stress for tau_j.')
            except:
                print('Failed to fetch data for {}'.format(elements))
                print('Make sure data file contains: ')
                print('1. "structure", \n'
                      '2. "Vn" or "a" or "b", and elastic moduli \n '
                      '3. "E_w_i"  \n ')
            self.data_of_ssmodels["BCC_screw_Suzuki_RWASM-2020"] = elements_data
        
