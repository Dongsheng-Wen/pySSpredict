import pandas as pd 
import copy
try:
    from make_prediction.read_input import  *
    from make_prediction.models import *
    from make_prediction.make_composition import build_mesh_ternary
except:
    from .read_input import *
    from .models import *
    from .make_composition import build_mesh_ternary



class ternary_calc_strength:

    def __init__(self,comp={},modelname=[]):

        self.comp = {}
        self.temp = np.array([300.])
        self.strain_r = 0.001
        self.conditions = None
        self.elements = []
        self.tcdatabase = None
        self.no_liquidus_T = True 
        self.T_l = None
        self.calc_data = []
        self.ssmodels = []
        self.calc_data_all = None
        self.all_elements = []
        self.ssmodels_all = ["FCC_Varvenne-Curtin-2016",
                            "BCC_edge_Maresca-Curtin-2019",
                            "BCC_screw_Maresca-Curtin-2019",
                            "BCC_screw_Suzuki_RWASM-2020"]

    def set_comp_mesh(self,from_dict=None,
        psA = None, psB = None, psC = None,
        comp_increment = None,
        alias=None,remove_T_l=True):
        # parameter: alloy composition
        # set the composition mesh for the ternary (or pseudo-ternary) system
        # supply a dictionary as follows and use from_dict
        # e.g. 
        ''' 
        {
            "increment": 1,
            "psA": {"grouped_elements":["Nb"],
                    "ratio": [1],
                    "range":[0,100]},
            "psB": {"grouped_elements":["Mo"],
                    "ratio": [1],
                    "range":[0,100]},
            "psC": {"grouped_elements":["W"],
                    "ratio": [1], 
                    "range":[0,100]} # 0 to 100 at.%
        }
        '''
        # if no from_dict is supply, psA, psB, psC, and comp_increment should be supplied 


        # alias = dict 
        # if more than one data set provided for one element
        # user can specify the alias of that atom 
        # e.g. {"Nb":"Nb1", "Ti":"Ti", "Zr":"Zr2"}

        if from_dict:
            self.comp_increment = from_dict['increment']
            self.psA_elements = from_dict['psA']['grouped_elements']
            self.psB_elements = from_dict['psB']['grouped_elements']
            self.psC_elements = from_dict['psC']['grouped_elements']
            self.psA_ratio = from_dict['psA']['ratio']
            self.psB_ratio = from_dict['psB']['ratio']
            self.psC_ratio = from_dict['psC']['ratio']
            self.psA_range = from_dict['psA']['range']
            self.psB_range = from_dict['psB']['range']
            self.psC_range = from_dict['psC']['range']

        else:
            self.comp_increment = comp_increment
            self.psA_elements = psA['grouped_elements']
            self.psB_elements = psB['grouped_elements']
            self.psC_elements = psC['grouped_elements']
            self.psA_ratio = psA['ratio']
            self.psB_ratio = psB['ratio']
            self.psC_ratio = psC['ratio']
            self.psA_range = psA['range']
            self.psB_range = psB['range']
            self.psC_range = psC['range']
        ps_tags = ['A','B','C']
        ps_elements = [self.psA_elements,self.psB_elements,self.psC_elements]
        ps_ranges = [self.psA_range,self.psB_range,self.psC_range]
        ps_ratios = [self.psA_ratio,self.psB_ratio,self.psC_ratio]
        self.mesh = build_mesh_ternary(self.comp_increment,
                          self.psA_range,self.psB_range,self.psC_range,
                          self.psA_elements,self.psB_elements,self.psC_elements,
                          self.psA_ratio,self.psB_ratio,self.psC_ratio)
        comp_tmp = self.mesh.make_mesh()
        
        self.comp_pst, self.comp_elements = self.mesh.comp_assign() #pd.dataframes
        
        self.comp = self.comp_elements/100 # pd.dataframes

        elements = self.comp_elements.columns.to_list()

        self.elements = elements
        self.alias = alias

        # comp_dict is for liquidus temperature calculations
        self.comp_dict = {}
        for e in self.elements:
            self.comp_dict[e] = np.array(self.comp[e])
            if e not in self.all_elements:
                self.all_elements.append(e)
        
        print_line = 'Alloy system: '
        for i in self.elements:
            print_line += '{} '.format(str(i))
        print(print_line)

        print('Ternary composition grid generated: ')
        for _vertex,_element,_range,_ratio in zip(ps_tags,ps_elements,ps_ranges,ps_ratios):
            if len(_element) > 1:
                _e_ratio = ''
                for _e,_ratio in zip(_element,_ratio):
                    _e_ratio = _e_ratio + '{}{}'.format(str(_e),str(_ratio))
                print('Vertex: {}; Vertex Elements: {}; Composition range: {}-{} at.%'.format(str(_vertex),str(_e_ratio),str(_range[0]),str(_range[1]) ) )
            else:
                print('Vertex: {}; Vertex Element: {}; Composition range: {}-{} at.%'.format(str(_vertex),str(_element[0]),str(_range[0]),str(_range[1]) ) )

        if remove_T_l:
            self.T_l = None

    def set_temp(self,temp):
        # temperature 1d array
        self.temp = np.atleast_1d(np.array(temp))
        self.conditions = [self.temp, self.strain_r]
    def set_strain_r(self,strain_r):
        # strain rate, float
        self.strain_r = np.float(strain_r)
        self.conditions = [self.temp, self.strain_r]
    def exp_conditions(self): 
        # set temperature and strain rate
        self.conditions = [self.temp, self.strain_r]
        print('Temperature: {}'.format(self.conditions[0]))
        print('Strain rate: {}'.format(self.conditions[1]))

    def set_T_l(self,T_l):
        self.T_l = np.atleast_1d(np.array(T_l)) 

    def grab_data(self,jsfh=None,from_dict=None):
        if len(self.ssmodels)!=0:
            # self.data_of_ssmodels: type = dict 
            # key = model name, value = dictionary containing elemental data for that model. 
            self.data_handle = get_elements_data(self.elements,self.ssmodels,fh=jsfh,from_dict=from_dict,alias=self.alias)
            self.data_of_ssmodels = self.data_handle.data_of_ssmodels
            
        else:
            print('No model selected.')
            print('Please set models first. Example: use_model("FCC_Varvenne-Curtin-2016"). ')


    def calc_T_liquidus(self,composition,tcdatabase=None):
        try:
            from make_prediction.tc_liquidus_wrapper import liquidus_wrapper
        except:
            from .tc_liquidus_wrapper import liquidus_wrapper
        print('Will calculate liquidus temperature using TC-Python.')
        self.change_tcdatabase(tcdatabase=tcdatabase)
        self.no_liquidus_T = False

        self.T_l_calculator = liquidus_wrapper(self.comp_dict, # dict of compositions
                                    tcdatabase=self.tcdatabase)
        self.T_l_calculator.calc_liquidus()
        self.T_l = (self.T_l_calculator.T_liquidus)
        print('Liquidus Temperature: {} K.'.format(self.T_l))

    def calc_phase_diagram(self,elements,temperature):
        try:
            from make_prediction.tc_isothermo_wrapper import ternary_isothermal_pd
        except:
            from .tc_isothermo_wrapper import ternary_isothermal_pd

        temperature = 1169
        self.ternary_pd_calc = ternary_isothermal_pd(elements,temperature)
        self.ternary_pd_calc.calc_diagram()
        self.ternary_pd_calc.make_diagram_data()
        self.pd_data = self.ternary_pd_calc.pd_data
        

    def set_adjustables(self,
                        f1=None,f2=None,alpha=None, # edge models
                        kink_width=None,Delta_V_p_para=None,Delta_E_p_para=None,# MC screw model
                        tau_i_exponent=None,dislocation_density=None,trial_kappa=None,trial_tau_k=None # Suzuki screw model
                        ):
        # make_adjustables: read_input.make_adjustables
        # set adjustable parameters for every model included.
        # self.adjustables: dict
        # key = model name
        # value = dict containing parameters
        adjustables_hd = make_adjustables(self.ssmodels,
                        f1=f1,f2=f2,alpha=alpha, # edge models
                        kink_width=kink_width,Delta_V_p_para=Delta_V_p_para,Delta_E_p_para=Delta_E_p_para,# MC screw model
                        tau_i_exponent=tau_i_exponent,dislocation_density=dislocation_density,trial_kappa=trial_kappa,trial_tau_k=trial_tau_k # Suzuki screw model
                        )
        self.adjustables = adjustables_hd.adjustables

    def calculate(self):
        if self.calc_data_all is not None:
            self.calc_data = [self.calc_data_all]
        if self.conditions is None:
            self.exp_conditions()
        
        for ssmodel in self.ssmodels:
            print('Preparing Calculation -> {}.'.format(ssmodel))
            elements_data = self.data_of_ssmodels[ssmodel]
            adjustables = self.adjustables[ssmodel]
            calc_done = []
            for i in range(len(self.comp)):
                compositions = pd.DataFrame(data=self.comp,index=[i]).reset_index(drop=True) * 100 # x100 because model functions take 0-100 at.%
                pst_comp = pd.DataFrame(data=self.comp_pst,index=[i]).reset_index(drop=True)
                
                if (ssmodel == 'FCC_Varvenne-Curtin-2016') and self.check_data(ssmodel=ssmodel): 
                    # single FCC edge model
                    structure = 'fcc'
                    
                    model_core = ss_edge_model_T(adjustables, # dislocation related parameters
                                        self.conditions, # exp_conditions = [[list of T],strain_r]
                                        compositions, # pd.df 
                                        elements_data, # dict 
                                        structure      # 'fcc' or 'bcc'
                                        ) 
                elif (ssmodel == 'BCC_edge_Maresca-Curtin-2019') and self.check_data(ssmodel=ssmodel): 
                    # single BCC edge model
                    structure = 'bcc'
                
                    model_core = ss_edge_model_T(adjustables, # dislocation related parameters
                                        self.conditions, # exp_conditions = [[list of T],strain_r]
                                        compositions, # pd.df 
                                        elements_data, # dict 
                                        structure      # 'fcc' or 'bcc'
                                        ) 

                elif (ssmodel == 'BCC_screw_Maresca-Curtin-2019') and self.check_data(ssmodel=ssmodel): 
                    # single BCC M-C screw model
                    structure = 'bcc'
                    model_core = ss_model_M_C_screw_single(adjustables, 
                                        self.conditions, # exp_conditions = [[list of T],strain_r]
                                        compositions,  # pd.df 
                                        elements_data   # dict 
                                        )

                elif (ssmodel == 'BCC_screw_Suzuki_RWASM-2020') and self.check_data(ssmodel=ssmodel): 
                    # single BCC Suzuki model
                    structure = 'bcc'
                    if self.T_l is not None:
                        T_liq = np.atleast_1d(np.array(self.T_l[i]))
                    else:
                        T_liq = None 

                    model_core = Suzuki_model_RWASM_jog_T(adjustables,
                                        self.conditions,
                                        compositions,
                                        elements_data,
                                        T_l = T_liq)
                else:
                    model_core = None

                if model_core is not None:
                    model_core.calculate()
                    res = model_core.calc_data_all
                    res['comp(psA)'] = self.comp_pst['comp(psA)'][i]
                    res['comp(psB)'] = self.comp_pst['comp(psB)'][i]
                    res['comp(psC)'] = self.comp_pst['comp(psC)'][i]
                    res['strain_rate'] = self.strain_r
                    res['model'] = ssmodel
                    res['adjustables'] = str(adjustables)
                    res['structure'] = structure
                    calc_done.append(True)
                    res['done'] = True
                    self.calc_data.append(res)
                else:
                    calc_done.append(False)
                    res = compositions
                    res['comp(psA)'] = self.comp_pst['comp(psA)'][i]
                    res['comp(psB)'] = self.comp_pst['comp(psB)'][i]
                    res['comp(psC)'] = self.comp_pst['comp(psC)'][i]
                    res['done'] = False
                    self.calc_data.append(res)
            if all(calc_done):
                print('Done Calculation. -> {}.'.format(ssmodel))
            else:
                print('WARNING: Some calculations maynot be performed. Check your input/output data -> {}.'.format(ssmodel))

        self.calc_data_all = pd.concat(self.calc_data,axis=0,ignore_index=True)
        try:
            self.calc_data_all = self.calc_data_all.drop_duplicates(subset=self.elements+['T','strain_rate','structure','jog_activated','model','adjustables'],keep='last')
        except:
            self.calc_data_all = self.calc_data_all.drop_duplicates()
        

    def change_tcdatabase(self,tcdatabase=None):
        
        if tcdatabase is not None:
            self.tcdatabase = tcdatabase
            #print('Change TC database: {}'.format(tcdatabase))
        elif tcdatabase is None and self.tcdatabase is not None:
            #print('TC database: {}'.format(tcdatabase))
            pass 
        else:
            self.tcdatabase = 'TCHEA4' 
            #print('No TC database specified. Use default TCHEA4.')

        print('Will use ThermoCalc database: {}'.format(self.tcdatabase))

    def use_model(self,modelname):
        # add model for prediction
        if modelname in self.ssmodels_all:
            if modelname not in self.ssmodels:
                self.ssmodels.append(modelname)
            else: 
                print('{} exists.'.format(modelname))
        else:
            print('Unknown model name. Supported models are: {}.'.format(self.ssmodels_all))
        print('Currently using model(s): {}'.format(self.ssmodels))


    def remove_model(self,modelname):
        try:
            self.ssmodels.remove(modelname)
            print('removing {} from the model set.'.format(modelname))
        except:
            print('{} not in the model list.'.format(modelname))
        print('Currently using model(s): {}.'.format(self.ssmodels))

    def reorganize_data(self,newcolumns=None):
        
        if newcolumns is None:
            data_cols = ['comp(psA)','comp(psB)','comp(psC)']
            for element_i in self.all_elements: 
                if element_i in self.calc_data_all.columns.to_list():
                    data_cols.append(element_i)
            if 'T_liquidus' in self.calc_data_all.columns.to_list():
                data_cols.append('T_liquidus')
                data_cols.append('jog_activated')
            pretty_data_cols = data_cols + ['T','strain_rate','tau_y','structure','model','adjustables']
        else: 
            pretty_data_cols = newcolumns
        
        self.pretty_calc_data = self.calc_data_all[pretty_data_cols]

    def check_data(self,ssmodel=None):
        data_ok = True
        if ssmodel is None:
            for ssmodel in self.ssmodels:
                data_check = self.data_of_ssmodels[ssmodel]
                if len(data_check) == len(self.elements):
                    for element in self.elements:
                        if not data_check[element]:
                            print('No data for {}. Make sure to have the relevant data for using the model -> {}.'.format(
                                element,ssmodel))
                            data_ok = False
                else: 
                    print('Not sufficient data for model -> {}'.format(ssmodel))
                    data_ok = False
        else: 
            data_check = self.data_of_ssmodels[ssmodel]
            if len(data_check) == len(self.elements):
                for element in self.elements:
                    if not data_check[element]:
                        data_ok = False
            else: 
                data_ok = False 
        return data_ok 