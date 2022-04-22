import pandas as pd 
import copy
try:
    from make_prediction.read_input import  *
    from make_prediction.models import *
except:
    from .read_input import *
    from .models import *



class single_calc_strength:

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

    def set_comp(self,comp,alias=None,remove_T_l=True):
        # parameter: alloy composition
        # comp = dict 
        # e.g. {"Nb":33, "Ti":33, "Zr":33}
        # alias = dict 
        # if more than one data set provided for one element
        # user can specify the alias of that atom 
        # e.g. {"Nb":"Nb1", "Ti":"Ti", "Zr":"Zr2"}
        try:
            #extract elements
            elements = [i for i in comp.keys()]
            #sum composition to ensure it is 100% or 1
            comp_sum = float(sum([float(comp[ele]) for ele in elements]))
            #if sum to 100, convert to mole fraction
            if comp_sum == 100.0:
                for ele in elements:
                    comp[ele] = float(comp[ele])/100.
            # check
            if (np.abs(comp_sum - 100.0) > 1E-8) and (np.abs(comp_sum - 1.0) > 1E-10):
                print('Your composition does not sum to 1.0. Perform normalization.')
                for ele in elements:
                    comp[ele] = float(comp[ele])/comp_sum
        except:
            print('Please provide a dictionary as input.')

        self.comp = comp
        self.elements = elements
        self.alias = alias
        for e in self.elements:
            if e not in self.all_elements:
                self.all_elements.append(e)

        print_line = 'Alloy system: '
        for i in comp.keys():
            print_line += '{}: {} '.format(str(i),str(round(comp[i],2)))
        print(print_line)
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


    def calc_T_liquidus(self,tcdatabase=None):
        try:
            from make_prediction.tc_liquidus_wrapper import liquidus_wrapper
        except:
            from .tc_liquidus_wrapper import liquidus_wrapper
        print('Will calculate liquidus temperature using TC-Python.')
        self.change_tcdatabase(tcdatabase=tcdatabase)
        self.no_liquidus_T = False
        self.T_l_calculator = liquidus_wrapper(self.comp, # dict 
                                    tcdatabase=self.tcdatabase)
        self.T_l_calculator.calc_liquidus()
        self.T_l = (self.T_l_calculator.T_liquidus)
        print('Liquidus Temperature: {} K.'.format(self.T_l))

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
        compositions = pd.DataFrame(data=self.comp,index=[0]) * 100 # x100 because model functions take 0-100 at.%
        if self.conditions is None:
            self.exp_conditions()
        
        for ssmodel in self.ssmodels:
            print('Preparing Calculation -> {}.'.format(ssmodel))
            elements_data = self.data_of_ssmodels[ssmodel]
            adjustables = self.adjustables[ssmodel]
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
                model_core = Suzuki_model_RWASM_jog_T(adjustables,
                                    self.conditions,
                                    compositions,
                                    elements_data,
                                    T_l = self.T_l)
            else:
                model_core = None

            if model_core is not None:
                model_core.calculate()
                print('Done Calculation. -> {}.'.format(ssmodel))
                res = model_core.calc_data_all
                res['strain_rate'] = self.strain_r
                res['model'] = ssmodel
                res['adjustables'] = str(adjustables)
                res['structure'] = structure
                self.calc_data.append(res)
            else:
                print('NO Calculation Performed. Check your data -> {}.'.format(ssmodel))
                self.calc_data.append(pd.DataFrame(columns=compositions.columns.to_list()))

        self.calc_data_all = pd.concat(self.calc_data,axis=0,ignore_index=True)
        try:
            self.calc_data_all = self.calc_data_all.drop_duplicates(subset=self.elements+['T','strain_rate','structure','jog_activated','model','adjustables'],keep='last')
        except:
            self.calc_data_all = self.calc_data_all.drop_duplicates()
        

    def change_tcdatabase(self,tcdatabase=None):
        
        if tcdatabase is None:
            self.tcdatabase = 'TCHEA4'
            print('No TC database specified. Use default TCHEA4.')
        else: 
            self.tcdatabase = tcdatabase

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
            data_cols = []
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