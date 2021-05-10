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
                self.mode = 'edge_single_calculation'
            else:
                print('Please indicate correct model and calculation type in the input file.\n'
                      'Use --format to check sample inputs.')
        elif self.model in model_BCC_screw_Curtin:
            if "pseudo-ternary" in self.data: 
                self.mode = 'screw_ternary'
            elif "properties" in self.data:
                self.mode = 'screw_single_calculation'
            else:
                print('Please indicate correct model and calculation type in the input file.\n'
                      'Use --format to check sample inputs.')
        elif self.model in model_BCC_screw_Suzuki: 
            if "pseudo-ternary" in self.data:
                self.mode = "screw_suzuki_ternary_calculation"
            elif "c" in self.data['elements'][list(self.data['elements'].keys())[0]]:
                self.mode = "screw_suzuki_single_calculation"
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
        def check_input_data(element_list,element_data):
            for element in element_list:
                if not element_data[element]:
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
        
        self.element_data = self.data['elements']
        
        self.elements_ABC = self.psA_elements + self.psB_elements + self.psC_elements
        
        self.temperature = self.data['conditions']['temperature']
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
        self.dislocation_properties = [self.alpha,self.f1,self.f2]
        for element_i in self.elements_ABC:
            # compute E, nu, G for elements if not supplied
            # two of the E/nu/G must be supplied to calculate the missing one
            if not 'nu' in self.element_data[element_i]: 
                self.element_data[element_i]['nu'] = round(self.element_data[element_i]['E']/2/self.element_data[element_i]['G'] - 1,3)
            if not 'G' in self.element_data[element_i]:
                self.element_data[element_i]['G'] = round(self.element_data[element_i]['E']/2/(self.element_data[element_i]['nu'] + 1),1)
            if not 'E' in self.element_data[element_i]:
                self.element_data[element_i]['E'] = round(self.element_data[element_i]['G']*2*(self.element_data[element_i]['nu'] + 1),1)
            
            # compute a/b/Vn for elements based on lattice structure
            # one of the a/b/Vn must be supplied to compute the other two
            if self.structure == 'fcc':
                # fcc: Vn = a^3/4, b = a/sqrt(2)
                if 'Vn' in self.element_data[element_i]:
                    self.element_data[element_i]['a'] = (self.element_data[element_i]['Vn']*4)**(1/3)
                    self.element_data[element_i]['b'] = self.element_data[element_i]['a']/np.sqrt(2)
                elif 'b' in self.element_data[element_i]:
                    self.element_data[element_i]['a'] = self.element_data[element_i]['b']*np.sqrt(2)
                    self.element_data[element_i]['Vn'] = self.element_data[element_i]['a']**3/4
                elif 'a' in self.element_data[element_i]:
                    self.element_data[element_i]['b'] = self.element_data[element_i]['a']/np.sqrt(2)
                    self.element_data[element_i]['Vn'] = self.element_data[element_i]['a']**3/4
            else: 
                # bcc: Vn = a^3/2, b = a*sqrt(3)/2
                if 'Vn' in self.element_data[element_i]:
                    self.element_data[element_i]['a'] = (self.element_data[element_i]['Vn']*2)**(1/3)
                    self.element_data[element_i]['b'] = self.element_data[element_i]['a']*np.sqrt(3)/2
                elif 'b' in self.element_data[element_i]:
                    self.element_data[element_i]['a'] = self.element_data[element_i]['b']*2/np.sqrt(3)
                    self.element_data[element_i]['Vn'] = self.element_data[element_i]['a']**3/2
                elif 'a' in self.element_data[element_i]:
                    self.element_data[element_i]['b'] = self.element_data[element_i]['a']*np.sqrt(3)/2
                    self.element_data[element_i]['Vn'] = self.element_data[element_i]['a']**3/2
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
        def check_input_data(element_list,element_data):
            for element in element_list:
                if not element_data[element]:
                    print('No data found: {}'.format(element))
                
    def grab_properties_curtin_edge(self):

        
        self.element_data = self.data['elements']
        
        self.elements_order = self.data['compositions']['element_order']
            
        self.concentrations = self.data['compositions']['concentrations']
        
        self.temperature_range = np.arange(self.data['conditions']['temperature']['min'],
                                    self.data['conditions']['temperature']['max']+self.data['conditions']['temperature']['inc'],
                                    self.data['conditions']['temperature']['inc'])
        element_composition = {}
        for i in range(len(self.elements_order)):
            c_x = np.array(self.concentrations).transpose()[i]
            c_x_T = (np.ones((len(self.temperature_range),len(c_x)))*np.array(c_x)).transpose().flatten()
            element_composition[self.elements_order[i]] = c_x_T
        self.temperature = (np.ones((len(self.concentrations),len(self.temperature_range))) * self.temperature_range).flatten()
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
        self.dislocation_properties = [self.alpha,self.f1,self.f2]
        for element_i in self.elements_order:
            # compute E, nu, G for elements if not supplied
            # two of the E/nu/G must be supplied to calculate the missing one
            if not 'nu' in self.element_data[element_i]: 
                self.element_data[element_i]['nu'] = round(self.element_data[element_i]['E']/2/self.element_data[element_i]['G'] - 1,3)
            if not 'G' in self.element_data[element_i]:
                self.element_data[element_i]['G'] = round(self.element_data[element_i]['E']/2/(self.element_data[element_i]['nu'] + 1),1)
            if not 'E' in self.element_data[element_i]:
                self.element_data[element_i]['E'] = round(self.element_data[element_i]['G']*2*(self.element_data[element_i]['nu'] + 1),1)
            
            # compute a/b/Vn for elements based on lattice structure
            # one of the a/b/Vn must be supplied to compute the other two
            if self.structure == 'fcc':
                # fcc: Vn = a^3/4, b = a/sqrt(2)
                if 'Vn' in self.element_data[element_i]:
                    self.element_data[element_i]['a'] = (self.element_data[element_i]['Vn']*4)**(1/3)
                    self.element_data[element_i]['b'] = self.element_data[element_i]['a']/np.sqrt(2)
                elif 'b' in self.element_data[element_i]:
                    self.element_data[element_i]['a'] = self.element_data[element_i]['b']*np.sqrt(2)
                    self.element_data[element_i]['Vn'] = self.element_data[element_i]['a']**3/4
                elif 'a' in self.element_data[element_i]:
                    self.element_data[element_i]['b'] = self.element_data[element_i]['a']/np.sqrt(2)
                    self.element_data[element_i]['Vn'] = self.element_data[element_i]['a']**3/4
            else: 
                # bcc: Vn = a^3/2, b = a*sqrt(3)/2
                if 'Vn' in self.element_data[element_i]:
                    self.element_data[element_i]['a'] = (self.element_data[element_i]['Vn']*2)**(1/3)
                    self.element_data[element_i]['b'] = self.element_data[element_i]['a']*np.sqrt(3)/2
                elif 'b' in self.element_data[element_i]:
                    self.element_data[element_i]['a'] = self.element_data[element_i]['b']*2/np.sqrt(3)
                    self.element_data[element_i]['Vn'] = self.element_data[element_i]['a']**3/2
                elif 'a' in self.element_data[element_i]:
                    self.element_data[element_i]['b'] = self.element_data[element_i]['a']*np.sqrt(3)/2
                    self.element_data[element_i]['Vn'] = self.element_data[element_i]['a']**3/2
                    
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
        self.element_data = self.data['elements']
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
        
        # adjustable scalers
        self.adjustable_scalers = self.data['adjustables']
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
        # adjustable scalers
        self.adjustable_scalers = self.data['adjustables']
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
        self.element_data = self.data['elements']
        # adjustable scalers
        self.adjustable_scalers = self.data['adjustables']
        # exp conditions
        self.experiment_conditions = self.data['conditions']
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
        self.element_data = self.data['elements']
        # adjustable scalers
        self.adjustable_scalers = self.data['adjustables']
        # exp conditions
        self.experiment_conditions = self.data['conditions']
        # output file
        try: 
            self.savefilename = self.data['savefile']
        except:
            self.savefilename = self.data["material"] + '_out'
        
        