from sspredict.make_prediction.tc_isothermo_wrapper import single_point_isothermo
from sspredict.make_prediction.single_calc_wrapper import single_calc_strength
from sspredict.make_prediction.make_composition import composition_generator
try:
    from tc_python import *
except:
    print('cannot import TCPython')
import os 
import concurrent.futures
import pandas as pd
import itertools 
import pandas as pd 
import numpy as np 


class parallel_wrapper:
    # perform TC python and solid solution stress
    # 1. for an arbitrary N-component alloy
    #    calculate the stable phases and solid solution stress
    def __init__(self,
                elements,   # elements you want to include 
                elemental_data, # dictionary of elemental data containing essential properties 
                N_component, # int: number of components, e.g. 4
                comp_inc, 
                temperature,
                tcdatabase,
                ss_model_parameters = {"strain_rate":0.001},
                solid_solution_model_name = ['BCC_edge_Maresca-Curtin-2019'], # a list of model(s)
                output_file_prefix = None):
        # tcdatabase=self.tcdatabase
        self.tcdatabase=tcdatabase
        # number of processors
        self.num_proc = os.cpu_count()
        
        #elements,   # elements you want to include 
        self.elements = elements
        
        #elemental_data, # dictionary of elemental data containing essential properties 
        self.elemental_data = elemental_data
        
        #solid_solution_model_name = ['BCC_edge_Maresca-Curtin-2019'] # a list of model(s)
        self.solid_solution_model_name = solid_solution_model_name
        
        #N_component, # int: number of components, e.g. 4
        self.N_component = N_component
        
        # comp_inc, # increment of elemental concentration in an alloy
        self.comp_inc = comp_inc
        
        # temperature
        self.temperature = temperature
        
        # ss_model_parameters: dict
        # must contain {"strain_rate":0.001}
        if ss_model_parameters is not None: 
            self.ss_model_parameters = ss_model_parameters
            
        self.output_file_prefix = 'stdout'
        if output_file_prefix is not None: 
            self.output_file_prefix = output_file_prefix
    def split(self, a, n):
        # divide a list into n parts with approximately the same length
        k, m = divmod(len(a), n)
        return [a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n)]
    
    def generate_alloy_list(self):
        comp_list = []
        alloy_systems = [a for a in itertools.combinations(self.elements,self.N_component)]
        self.comp_gen = composition_generator(self.comp_inc,self.N_component)
        self.comp_gen.get_all_comps()
        # split the list into num_proc slices
        self.alloy_systems_slices = self.split(alloy_systems,self.num_proc)
        print('Total Number of Alloy Systems: {} Divided Into {} Groups'.format(len(alloy_systems),self.num_proc))
    
    def calc_slice(self,slice_index):
        single_calc = single_calc_strength()
        for model in self.solid_solution_model_name:
            single_calc.use_model(model)
        # setup parameters
        f1=None;f2=None;alpha=None, # edge models
        kink_width=None;Delta_V_p_para=None;Delta_E_p_para=None # MC screw model
        tau_i_exponent=None;dislocation_density=None;trial_kappa=None;trial_tau_k=None # Suzuki screw model
        params = [None,None,None,None,None,None,None,None,None,None,None]
        param_names = ['f1','f2','alpha',
                       'kink_width','Delta_V_p_para','Delta_E_p_para',
                       'tau_i_exponent','dislocation_density','trial_kappa','trial_tau_k']
        
        for param in self.ss_model_parameters.keys():
            try: 
                param_index = param_names.index(param)
                params[param_index] = self.ss_model_parameters[key]
            except:
                pass
        
        single_calc.set_adjustables(f1=params[0],f2=params[1],alpha=params[2], # edge models
        kink_width=params[3],Delta_V_p_para=params[4],Delta_E_p_para=params[5], # MC screw model
        tau_i_exponent=params[6],dislocation_density=params[7],trial_kappa=params[8],trial_tau_k=params[9]) # Suzuki screw model
            
        single_calc.set_temp(self.temperature)
        single_calc.set_strain_r(self.ss_model_parameters["strain_rate"])
        single_calc.exp_conditions()
        tc_single_calc = single_point_isothermo(temperature=self.temperature,tcdatabase=self.tcdatabase)
        
        for alloy_sys in self.alloy_systems_slices[slice_index]:
            eles = (dict(zip(alloy_sys,np.ones(len(alloy_sys)))))
            with TCPython() as Setup_TC:
                tc_single_calc.set_system(Setup_TC=Setup_TC,comp=eles)    
                for comp_i in self.comp_gen.comps:
                    comp = (dict(zip(alloy_sys,np.array(comp_i))))
                    print(comp)
                    # setup TC calculator 
                    # calculate equilibrium phases
                    tc_single_calc.set_comp(comp=comp)
                    tc_single_calc.calculate()
                    print(tc_single_calc.stable_phases)
                    # calculate strength
                    single_calc.set_comp(comp)
                    single_calc.grab_data(from_dict=elements_dict)
                    single_calc.calculate(stable_phases=str(tc_single_calc.stable_phases))
                    #single_calc.calc_data_all.to_csv('./{}_{}.csv'.format(self.output_file_prefix,str(int(slice_index))),sep=',')
        return single_calc.calc_data_all
    def run(self):
        self.parallel_data_list = []
        
        with concurrent.futures.ProcessPoolExecutor(self.num_proc) as executor:
            for result_from_process in zip(range(self.num_proc), executor.map(self.calc_slice, range(self.num_proc))):
                # params can be used to identify the process and its parameters
                parallel_data = result_from_process
                self.parallel_data_list.append(parallel_data)
            self.parallel_data_all = pd.concat(self.parallel_data_list,axis=0,ignore_index=True)
            self.parallel_data_all.to_csv('./{}_all.csv'.format(self.output_file_prefix),sep=',')
        print("Done")