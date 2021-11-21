import scipy 
import numpy as np
import pandas as pd
from scipy import stats 
import copy
import scipy.optimize as optimize
import scipy.integrate as integrate

class ss_edge_model_T:
# calculate solid solution strengthening contribution for FCC/BCC CCAs
# Edge dislocation models
# FCC model: Varvenne-Leyson-Ghazisaeidi-Curtin 2016: http://dx.doi.org/10.1016/j.actamat.2016.09.046
# BCC model: Maresca-Curtin 2020: https://doi.org/10.1016/j.actamat.2019.10.015
# for simeple calculations 
    def __init__(self,
                 adjustable_scalers,
                 exp_conditions,
                 comp_elements,
                 elements_data,
                 structure):
        
        # dislocation properties, alpha, f1, and f2
        try:
            self.alpha = float(adjustable_scalers[0])
            self.f_tau = float(adjustable_scalers[1])
            self.f_dEb = float(adjustable_scalers[2])
        except:
            self.alpha = float(adjustable_scalers['alpha'])
            self.f_tau = float(adjustable_scalers['f1'])
            self.f_dEb = float(adjustable_scalers['f2'])
        
        # experiment conditions, T, strain rate
        self.T_range = np.array(exp_conditions[0])
        self.ep = float(exp_conditions[1])
        self.ep0 = 10**4 #reference strain rate (/s)
        
        # some constants 
        self.boltzmann_J = 1.38064852*10**(-23) #J/K
        self.J2eV=6.2415093433*10**18 # covert J to eV
        
        # elemental data
        self.elements_order = comp_elements.columns.tolist()
        self.compositions_all = comp_elements #pandas df
        self.elements_data = copy.deepcopy(elements_data) #json
        
        # convert unit for properties
        # Vn: Å^3 to m^3
        # b: Å to m
        # a: Å to m
        # E: GPa to Pa
        # G: GPa to Pa
        
        for element_i in self.elements_order:
            self.elements_data[element_i]['Vn'] = elements_data[element_i]['Vn']*10**(-30)
            self.elements_data[element_i]['b'] = elements_data[element_i]['b']*10**(-10)
            self.elements_data[element_i]['a'] = elements_data[element_i]['a']*10**(-10)
            self.elements_data[element_i]['E'] = elements_data[element_i]['E']*10**(9)
            self.elements_data[element_i]['G'] = elements_data[element_i]['G']*10**(9)

        
        
        self.structure = structure
    
    def change_elements_data(self,elements_data):
        self.elements_data = elements_data
        for element_i in self.elements_order:
            self.elements_data[element_i]['Vn'] = elements_data[element_i]['Vn']*10**(-30)
            self.elements_data[element_i]['b'] = elements_data[element_i]['b']*10**(-10)
            self.elements_data[element_i]['a'] = elements_data[element_i]['a']*10**(-10)
            self.elements_data[element_i]['E'] = elements_data[element_i]['E']*10**(9)
            self.elements_data[element_i]['G'] = elements_data[element_i]['G']*10**(9)


    def FCC_V_L_G_C_2016_analytical(self):
        # FCC model: Varvenne-Leyson-Ghazisaeidi-Curtin 2016: http://dx.doi.org/10.1016/j.actamat.2016.09.046
        
        self.prefac_ty0 = 0.051
        self.Taylor_fac = 3.06
        self.prefac_dEb = 0.274
        # averaged properties
        cn_Vn = []
        cn_nu = []
        cn_G = []
        cn_b = []
        cn_E = []
        for element_i in self.elements_order:
            cn_Vn.append(self.compositions[element_i]/100*self.elements_data[element_i]['Vn'])
            cn_nu.append(self.compositions[element_i]/100*self.elements_data[element_i]['nu'])
            cn_G.append(self.compositions[element_i]/100*self.elements_data[element_i]['G'])
            cn_b.append(self.compositions[element_i]/100*self.elements_data[element_i]['b'])
            cn_E.append(self.compositions[element_i]/100*self.elements_data[element_i]['E'])

        self.aver_E = sum(cn_E);
        self.aver_V = sum(cn_Vn);
        self.aver_G = sum(cn_G)
        self.aver_Nu = sum(cn_nu)
        self.aver_b = sum(cn_b)
        
        i = 0;cn_Delta_Vn2=[]
        for element_i in self.elements_order:
            cn_Delta_Vn2.append(self.compositions[element_i]/100*(self.elements_data[element_i]['Vn']-self.aver_V)**2)
            
        self.sum_cndVn_b6 = sum(cn_Delta_Vn2)/self.aver_b**6;
        q_nu = ((1 + self.aver_Nu)/(1 - self.aver_Nu))
        
        self.dEb = self.prefac_dEb * self.f_dEb * self.alpha**(1/3)  * self.aver_G * self.aver_b**3   * q_nu**(2/3) * self.sum_cndVn_b6**(1/3)
        self.Ty0 = self.prefac_ty0 * self.f_tau * self.alpha**(-1/3) * self.aver_G * q_nu**(4/3) * self.sum_cndVn_b6**(2/3)
        self.Ty0_pc = self.Taylor_fac * self.Ty0
        delta_ss_low_T = self.Ty0 * (1 - ((self.boltzmann_J*self.T)/(self.dEb) * np.log(self.ep0/self.ep))**(2/3) )
        delta_ss_high_T = self.Ty0 * np.exp(-1/0.55 * self.boltzmann_J*self.T/self.dEb * np.log(self.ep0/self.ep) )
        self.delta_ss_low_T = delta_ss_low_T
        self.delta_ss_high_T = delta_ss_high_T
        Ty_threshold = self.Ty0/2
        self.delta_ss = self.Taylor_fac*np.array(delta_ss_low_T if delta_ss_low_T>=Ty_threshold else delta_ss_high_T)
        
        
    def BCC_M_C_2020_analytical(self):
        # BCC model: Maresca-Curtin-2019: https://doi.org/10.1016/j.actamat.2019.10.015
        
        self.prefac_ty0 = 0.051
        self.Taylor_fac = 3.06
        self.prefac_dEb = 0.274
        # averaged properties
        cn_Vn = []
        cn_nu = []
        cn_G = []
        cn_b = []
        cn_E = []
        for element_i in self.elements_order:
            cn_Vn.append(self.compositions[element_i]/100*self.elements_data[element_i]['Vn'])
            cn_nu.append(self.compositions[element_i]/100*self.elements_data[element_i]['nu'])
            cn_G.append(self.compositions[element_i]/100*self.elements_data[element_i]['G'])
            cn_b.append(self.compositions[element_i]/100*self.elements_data[element_i]['b'])
            cn_E.append(self.compositions[element_i]/100*self.elements_data[element_i]['E'])

        self.aver_E = sum(cn_E);
        self.aver_V = sum(cn_Vn);
        self.aver_G = sum(cn_G)
        self.aver_Nu = sum(cn_nu)
        self.aver_b = sum(cn_b)
        
        i = 0;cn_Delta_Vn2=[]
        for element_i in self.elements_order:
            cn_Delta_Vn2.append(self.compositions[element_i]/100*(self.elements_data[element_i]['Vn']-self.aver_V)**2)
            
        self.sum_cndVn_b6 = sum(cn_Delta_Vn2)/self.aver_b**6;
        q_nu = ((1 + self.aver_Nu)/(1 - self.aver_Nu))
        
        self.dEb = self.prefac_dEb * self.f_dEb * self.alpha**(1/3)  * self.aver_G * self.aver_b**3   * q_nu**(2/3) * self.sum_cndVn_b6**(1/3)
        self.Ty0 = self.prefac_ty0 * self.f_tau * self.alpha**(-1/3) * self.aver_G * q_nu**(4/3) * self.sum_cndVn_b6**(2/3)
        self.Ty0_pc = self.Taylor_fac * self.Ty0
        delta_ss_low_T = self.Ty0 * (1 - ((self.boltzmann_J*self.T)/(self.dEb) * np.log(self.ep0/self.ep))**(2/3) )
        delta_ss_high_T = self.Ty0 * np.exp(-1/0.55 * self.boltzmann_J*self.T/self.dEb * np.log(self.ep0/self.ep) )
        Ty_threshold = self.Ty0/2
        self.delta_ss_low_T = delta_ss_low_T
        self.delta_ss_high_T = delta_ss_high_T
        self.delta_ss = self.Taylor_fac*np.array(delta_ss_low_T if delta_ss_low_T>=Ty_threshold else delta_ss_high_T)
        
    def calculate(self):
        calc_header = ['T','V_ave','b_ave','E_ave','G_ave','nu_ave','sum_cnVn^2_b6','Ty0','Delta_Eb','tau_y','Delta_sigma_ss','tau_y_lT','tau_y_hT']
        element_header = self.compositions_all.columns.to_list()
        calc_data_all = pd.DataFrame(columns=element_header+calc_header)
        for composition_idx in range(len(self.compositions_all)):
            self.compositions = self.compositions_all.iloc[composition_idx]
           # print(self.compositions)
            for T in self.T_range:
                self.T = T
                #print(T)
                if self.structure == 'fcc':
                    self.FCC_V_L_G_C_2016_analytical()
                elif self.structure == 'bcc':
                    self.BCC_M_C_2020_analytical()
                
                calc_values = [self.T,
                               np.round(self.aver_V*10**30,2),
                               np.round(self.aver_b*10**10,4),
                               np.round( self.aver_E/10**9,2),
                               np.round( self.aver_G/10**9,2),
                               np.round( self.aver_Nu,4),
                               np.round(self.sum_cndVn_b6,8),
                               np.round(self.Ty0/10**6,2),
                               np.round(self.dEb*self.J2eV,4),
                               np.round(self.delta_ss/10**6/self.Taylor_fac,2),
                               np.round(self.delta_ss/10**6,2),
                               np.round(self.delta_ss_low_T/10**6,2),
                               np.round(self.delta_ss_high_T/10**6,2)
                               ]
                #print(calc_values)
                to_append = {**dict(zip(element_header,self.compositions.values)),**dict(zip(calc_header,calc_values))}
                #print(to_append)
                calc_data_all = calc_data_all.append(to_append,ignore_index=True)
            #print(to_append)
        self.calc_data_all = (calc_data_all)
        #print(self.calc_data_all)
    def writedata(self,foutname):
        outputfile = foutname 
        self.calc_data_all.to_csv(outputfile,index=False,sep='\t',float_format='%.2f')
        print('-------------------------------Data Saved to {}-------------------------------'.format(outputfile))





class ss_edge_model_batch_w_uncertainty:
# calculate solid solution strengthening contribution for FCC/BCC CCAs
# different from ss_edge_model, 
# consider the uncertainties in the elemental data input, lattice constants and elastic constants
# Edge dislocation models
# FCC model: Varvenne-Leyson-Ghazisaeidi-Curtin 2016: http://dx.doi.org/10.1016/j.actamat.2016.09.046
# BCC model: Maresca-Curtin 2019: https://doi.org/10.1016/j.actamat.2019.10.015
    def __init__(self,
                 ss_edge_model_batch,
                 dislocation_properties,
                 exp_conditions,
                 comp_elements,
                 comp_pst,
                 elements_data,
                 uncertainty_levels,
                 structure):
        

        self.dislocation_properties = dislocation_properties
        self.exp_conditions = exp_conditions
        self.T_range = np.array(self.exp_conditions[0])
        self.comp_elements = comp_elements
        self.elements_order = comp_elements.columns.tolist()
        self.comp_pst = comp_pst
        self.elements_data_save = copy.deepcopy(elements_data)
        self.structure = structure 
        self.J2eV=6.2415093433*10**18 # covert J to eV


        # uncertainty_levels controls the distribution of random variables of inputs
        # uncertainty_levels[0] for lattice constant a
        # uncertainty_levels[1] for elastic constants E, G, nu
        # use normal distribution
        # so uncertainty level is converted to standard deviation 
        # then the uncertainty will propagate to the predicted quantities. 
        # predicted quantity uncertainty will be appended to the predicted data. 
        
        self.uncertainty_levels = uncertainty_levels 
        
    def gen_rv(self):
        new_elements_data = copy.deepcopy(self.elements_data_save)
        if self.structure == 'fcc':
            for element_i in self.elements_order:
                new_elements_data[element_i]['a']  = stats.norm.rvs( self.elements_data_save[element_i]['a'],scale=self.elements_data_save[element_i]['a']*self.uncertainty_levels[0])
                new_elements_data[element_i]['b']  = new_elements_data[element_i]['b']
                new_elements_data[element_i]['Vn'] = new_elements_data[element_i]['Vn']
                new_elements_data[element_i]['E']  = stats.norm.rvs( self.elements_data_save[element_i]['E'],scale=self.elements_data_save[element_i]['E']*self.uncertainty_levels[1])
                new_elements_data[element_i]['G']  = stats.norm.rvs( self.elements_data_save[element_i]['G'],scale=self.elements_data_save[element_i]['G']*self.uncertainty_levels[1])
                new_elements_data[element_i]['nu']  = stats.norm.rvs( self.elements_data_save[element_i]['nu'],scale=self.elements_data_save[element_i]['nu']*self.uncertainty_levels[1])
        else:
            for element_i in self.elements_order:
                new_elements_data[element_i]['a']  = stats.norm.rvs( self.elements_data_save[element_i]['a'],scale=self.elements_data_save[element_i]['a']*self.uncertainty_levels[0])
                new_elements_data[element_i]['b']  = new_elements_data[element_i]['b']
                new_elements_data[element_i]['Vn'] = new_elements_data[element_i]['Vn']
                new_elements_data[element_i]['E']  = stats.norm.rvs( self.elements_data_save[element_i]['E'],scale=self.elements_data_save[element_i]['E']*self.uncertainty_levels[1])
                new_elements_data[element_i]['G']  = stats.norm.rvs( self.elements_data_save[element_i]['G'],scale=self.elements_data_save[element_i]['G']*self.uncertainty_levels[1])
                new_elements_data[element_i]['nu']  = stats.norm.rvs( self.elements_data_save[element_i]['nu'],scale=self.elements_data_save[element_i]['nu']*self.uncertainty_levels[1])
        return new_elements_data

    def calculate(self):
        self.calc_data_all = {}
        for T in self.T_range:
            self.T = T 
            self.exp_condition = []
            self.exp_condition.append([T])
            self.exp_condition.append(self.exp_conditions[1])
            self.aver_V_list = []
            self.aver_b_list = []
            self.aver_E_list = []
            self.aver_G_list = []
            self.aver_Nu_list = []
            self.sum_cndVn_b6_list = []
            self.Ty0_list = []
            self.dEb_list = []
            self.delta_ss_list = []
            i=0
            while i <=1000:

                self.elements_data = self.gen_rv()
                self.model = ss_edge_model_batch(self.dislocation_properties,
                                           self.exp_condition,
                                           self.comp_elements,self.comp_pst,
                                           self.elements_data,
                                           self.structure
                                          )
                #calculate data

                if self.structure == 'fcc':
                    self.model.FCC_V_L_G_C_2016_analytical()
                elif self.structure == 'bcc':
                    self.model.BCC_M_C_2020_analytical()

                self.aver_V_list.append(self.model.aver_V)
                self.aver_b_list.append(self.model.aver_b)
                self.aver_E_list.append(self.model.aver_E)
                self.aver_G_list.append(self.model.aver_G)
                self.aver_Nu_list.append(self.model.aver_Nu)
                self.sum_cndVn_b6_list.append(self.model.sum_cndVn_b6)
                self.Ty0_list.append(self.model.Ty0)
                self.dEb_list.append(self.model.dEb)
                self.delta_ss_list.append(self.model.delta_ss)
                i+=1


            self.aver_V = np.mean( np.array([ aver_V for aver_V in self.aver_V_list ]), axis=0 )
            self.aver_b = np.mean( np.array([ aver_b for aver_b in self.aver_b_list ]), axis=0 )
            self.aver_E = np.mean( np.array([ aver_E for aver_E in self.aver_E_list ]), axis=0 )
            self.aver_G = np.mean( np.array([ aver_G for aver_G in self.aver_G_list ]), axis=0 )
            self.aver_Nu = np.mean( np.array([ aver_Nu for aver_Nu in self.aver_Nu_list ]), axis=0 )
            self.aver_sum_cndVn_b6 = np.mean( np.array([ sum_cndVn_b6 for sum_cndVn_b6 in self.sum_cndVn_b6_list ]), axis=0 )
            self.aver_Ty0 = np.mean( np.array([ Ty0 for Ty0 in self.Ty0_list ]), axis=0 )
            self.aver_dEb = np.mean( np.array([ dEb for dEb in self.dEb_list ]), axis=0 )
            self.aver_delta_ss = np.mean( np.array([delta_ss for delta_ss in self.delta_ss_list ]), axis=0 )
            # evaluate uncertainty standard deviation
            self.std_V = np.std( np.array([ aver_V for aver_V in self.aver_V_list ]), axis=0 )
            self.std_b = np.std( np.array([ aver_b for aver_b in self.aver_b_list ]), axis=0 )
            self.std_E = np.std( np.array([ aver_E for aver_E in self.aver_E_list ]), axis=0 )
            self.std_G = np.std( np.array([ aver_G for aver_G in self.aver_G_list ]), axis=0 )
            self.std_Nu = np.std( np.array([ aver_Nu for aver_Nu in self.aver_Nu_list ]), axis=0 )
            self.std_sum_cndVn_b6 = np.std( np.array([ sum_cndVn_b6 for sum_cndVn_b6 in self.sum_cndVn_b6_list ]), axis=0 )
            self.std_Ty0 = np.std( np.array([ Ty0 for Ty0 in self.Ty0_list ]), axis=0 )
            self.std_dEb = np.std( np.array([ dEb for dEb in self.dEb_list ]), axis=0 )
            self.std_delta_ss = np.std( np.array([ delta_ss for delta_ss in self.delta_ss_list ]), axis=0 )


            self.calc_data = copy.deepcopy(comp_pst)
            self.calc_data['V_ave'] = self.aver_V*10**30
            self.calc_data['b_ave'] = self.aver_b*10**10
            self.calc_data['E_ave'] = self.aver_E/10**9
            self.calc_data['G_ave'] = self.aver_G/10**9
            self.calc_data['nu_ave'] = self.aver_Nu
            self.calc_data['T'] = [self.T for i in range(len(self.aver_Ty0))]
            self.calc_data['sum_cnVn^2_b6'] = [round(self.aver_sum_cndVn_b6[i],8) for i in range(len(self.aver_Ty0))]
            self.calc_data['Ty0'] = [round(self.aver_Ty0[i]/10**6,2) for i in range(len(self.aver_Ty0))]
            self.calc_data['Delta_Eb'] = [round(self.aver_dEb[i]*self.J2eV,4) for i in range(len(self.aver_Ty0))]
            self.calc_data['Delta_sigma_ss'] = [round(self.aver_delta_ss[i]/10**6,2) for i in range(len(self.aver_Ty0))]



            self.calc_data['std_V_ave'] = self.std_V*10**30
            self.calc_data['std_b_ave'] = self.std_b*10**10
            self.calc_data['std_E_ave'] = self.std_E/10**9
            self.calc_data['std_G_ave'] = self.std_G/10**9
            self.calc_data['std_nu_ave'] = self.std_Nu
            self.calc_data['std_sum_cnVn^2_b6'] = [round(self.std_sum_cndVn_b6[i],8) for i in range(len(self.aver_Ty0))] 
            self.calc_data['std_Ty0'] =  [round(self.std_Ty0[i]/10**6,2) for i in range(len(self.aver_Ty0))] 
            self.calc_data['std_Delta_Eb'] = [round(self.std_dEb[i]*self.J2eV,4) for i in range(len(self.aver_Ty0))] 
            self.calc_data['std_Delta_sigma_ss'] = [round(self.std_delta_ss[i]/10**6,2) for i in range(len(self.aver_Ty0))] 
            self.calc_data_all[str(self.T)] = self.calc_data
            
    def writedata(self,foutname):
        for key_i in self.calc_data_all.keys():
            calc_data = self.calc_data_all[key_i]
            outputfile = foutname+'_'+str(key_i)
            calc_data.to_csv(outputfile,index=False,sep='\t',float_format='%.2f')
            print('-------------------------------Data Saved to {}-------------------------------'.format(outputfile))






class ss_edge_model_T_w_uncertainty:
# calculate solid solution strengthening contribution for FCC/BCC CCAs
# different from ss_edge_model, 
# consider the uncertainties in the elemental data input, lattice constants and elastic constants
# Edge dislocation models
# FCC model: Varvenne-Leyson-Ghazisaeidi-Curtin 2016: http://dx.doi.org/10.1016/j.actamat.2016.09.046
# BCC model: Maresca-Curtin 2020: 
    def __init__(self,
                 ss_edge_model_T,
                 adjustable_scalers,
                 exp_conditions,
                 comp_elements,
                 elements_data,
                 uncertainty_levels,
                 structure):
        

        self.adjustable_scalers = adjustable_scalers
        self.exp_conditions = exp_conditions
        self.T_range = np.array(exp_conditions[0])
        self.compositions_all = comp_elements
        self.elements_order = comp_elements.columns.tolist()
        self.elements_data_save = copy.deepcopy(elements_data)
        self.structure = structure 
        self.J2eV=6.2415093433*10**18 # covert J to eV
        self.elements_data = elements_data

        # uncertainty_levels controls the distribution of random variables of inputs
        # uncertainty_levels[0] for lattice constant a
        # uncertainty_levels[1] for elastic constants E, G, nu
        # use normal distribution
        # so uncertainty level is converted to standard deviation 
        # then the uncertainty will propagate to the predicted quantities. 
        # predicted quantity uncertainty will be appended to the predicted data. 
        
        self.uncertainty_levels = uncertainty_levels 
        
    def gen_rv(self):

        new_elements_data = copy.deepcopy(self.elements_data_save)
        if self.structure == 'fcc':
            for element_i in self.elements_order:
                new_elements_data[element_i]['a']  = stats.norm.rvs( self.elements_data_save[element_i]['a'],scale=self.elements_data_save[element_i]['a']*self.uncertainty_levels[0])
                new_elements_data[element_i]['b']  = new_elements_data[element_i]['b']
                new_elements_data[element_i]['Vn'] = new_elements_data[element_i]['Vn']
                new_elements_data[element_i]['E']  = stats.norm.rvs( self.elements_data_save[element_i]['E'],scale=self.elements_data_save[element_i]['E']*self.uncertainty_levels[1])
                new_elements_data[element_i]['G']  = stats.norm.rvs( self.elements_data_save[element_i]['G'],scale=self.elements_data_save[element_i]['G']*self.uncertainty_levels[1])
                new_elements_data[element_i]['nu']  = stats.norm.rvs( self.elements_data_save[element_i]['nu'],scale=self.elements_data_save[element_i]['nu']*self.uncertainty_levels[1])
        else:
            for element_i in self.elements_order:
                new_elements_data[element_i]['a']  = stats.norm.rvs( self.elements_data_save[element_i]['a'],scale=self.elements_data_save[element_i]['a']*self.uncertainty_levels[0])
                new_elements_data[element_i]['b']  = new_elements_data[element_i]['b']
                new_elements_data[element_i]['Vn'] = new_elements_data[element_i]['Vn']
                new_elements_data[element_i]['E']  = stats.norm.rvs( self.elements_data_save[element_i]['E'],scale=self.elements_data_save[element_i]['E']*self.uncertainty_levels[1])
                new_elements_data[element_i]['G']  = stats.norm.rvs( self.elements_data_save[element_i]['G'],scale=self.elements_data_save[element_i]['G']*self.uncertainty_levels[1])
                new_elements_data[element_i]['nu']  = stats.norm.rvs( self.elements_data_save[element_i]['nu'],scale=self.elements_data_save[element_i]['nu']*self.uncertainty_levels[1])
        return new_elements_data

    def calculate(self):
        

        self.model = ss_edge_model_T(self.adjustable_scalers,
                           self.exp_conditions,
                           self.compositions_all,
                           self.elements_data,
                           self.structure
                          )
        calc_header = ['T','V_ave','b_ave','E_ave','G_ave','nu_ave','sum_cnVn^2_b6','Ty0','Delta_Eb','tau_y','Delta_sigma_ss',
                       'std_V_ave','std_b_ave','std_E_ave','std_G_ave','std_nu_ave','std_sum_cnVn','std_Ty0','std_Delta_Eb',
                      'std_tau_y','std_Delta_sigma_ss']
        element_header = self.compositions_all.columns.to_list()
        calc_data_all = pd.DataFrame(columns=element_header+calc_header)
        for composition_idx in range(len(self.compositions_all)):
            self.model.compositions = self.compositions_all.iloc[composition_idx]
               # print(self.compositions)
            for T in self.T_range:
                self.model.T = T
                self.aver_V_list = []
                self.aver_b_list = []
                self.aver_E_list = []
                self.aver_G_list = []
                self.aver_Nu_list = []
                self.sum_cndVn_b6_list = []
                self.Ty0_list = []
                self.dEb_list = []
                self.delta_ss_list = []
                i=0
                self.model.change_elements_data(copy.deepcopy(self.elements_data_save))
                #print(self.model.elements_data)
                if self.structure == 'fcc':
                    self.model.FCC_V_L_G_C_2016_analytical()
                elif self.structure == 'bcc':
                    self.model.BCC_M_C_2020_analytical()
                        
                self.aver_V = self.model.aver_V 
                self.aver_b = self.model.aver_b
                self.aver_E = self.model.aver_E 
                self.aver_G = self.model.aver_G 
                self.aver_Nu = self.model.aver_Nu
                self.aver_sum_cndVn_b6 = self.model.sum_cndVn_b6 
                self.aver_Ty0 = self.model.Ty0
                self.aver_dEb = self.model.dEb
                #self.aver_tau_y = self.model.tau_y 
                self.aver_delta_ss = self.model.delta_ss 
                
                while i <=200:

                    #self.elements_data = self.gen_rv()
                    self.model.change_elements_data(self.gen_rv())
                    if self.structure == 'fcc':
                        self.model.FCC_V_L_G_C_2016_analytical()
                    elif self.structure == 'bcc':
                        self.model.BCC_M_C_2020_analytical()

                    self.aver_V_list.append(self.model.aver_V)
                    self.aver_b_list.append(self.model.aver_b)
                    self.aver_E_list.append(self.model.aver_E)
                    self.aver_G_list.append(self.model.aver_G)
                    self.aver_Nu_list.append(self.model.aver_Nu)
                    self.sum_cndVn_b6_list.append(self.model.sum_cndVn_b6)
                    self.Ty0_list.append(self.model.Ty0)
                    self.dEb_list.append(self.model.dEb)
                    self.delta_ss_list.append(self.model.delta_ss)
                    i+=1

                ''' or you can take the average 
                self.aver_V = np.mean( np.array([ aver_V for aver_V in self.aver_V_list ]), axis=0 )
                self.aver_b = np.mean( np.array([ aver_b for aver_b in self.aver_b_list ]), axis=0 )
                self.aver_E = np.mean( np.array([ aver_E for aver_E in self.aver_E_list ]), axis=0 )
                self.aver_G = np.mean( np.array([ aver_G for aver_G in self.aver_G_list ]), axis=0 )
                self.aver_Nu = np.mean( np.array([ aver_Nu for aver_Nu in self.aver_Nu_list ]), axis=0 )
                self.aver_sum_cndVn_b6 = np.mean( np.array([ sum_cndVn_b6 for sum_cndVn_b6 in self.sum_cndVn_b6_list ]), axis=0 )
                self.aver_Ty0 = np.mean( np.array([ Ty0 for Ty0 in self.Ty0_list ]), axis=0 )
                self.aver_dEb = np.mean( np.array([ dEb for dEb in self.dEb_list ]), axis=0 )
                self.aver_delta_ss = np.mean( np.array([delta_ss for delta_ss in self.delta_ss_list ]), axis=0 )'''

                self.std_V = np.std( np.array([ aver_V for aver_V in self.aver_V_list ]), axis=0 )
                self.std_b = np.std( np.array([ aver_b for aver_b in self.aver_b_list ]), axis=0 )
                self.std_E = np.std( np.array([ aver_E for aver_E in self.aver_E_list ]), axis=0 )
                self.std_G = np.std( np.array([ aver_G for aver_G in self.aver_G_list ]), axis=0 )
                self.std_Nu = np.std( np.array([ aver_Nu for aver_Nu in self.aver_Nu_list ]), axis=0 )
                self.std_sum_cndVn_b6 = np.std( np.array([ sum_cndVn_b6 for sum_cndVn_b6 in self.sum_cndVn_b6_list ]), axis=0 )
                self.std_Ty0 = np.std( np.array([ Ty0 for Ty0 in self.Ty0_list ]), axis=0 )
                self.std_dEb = np.std( np.array([ dEb for dEb in self.dEb_list ]), axis=0 )
                self.std_delta_ss = np.std( np.array([ delta_ss for delta_ss in self.delta_ss_list ]), axis=0 )
                
                
                calc_values = [T,
                               np.round(self.aver_V*10**30,2),
                               np.round(self.aver_b*10**10,4),
                               np.round( self.aver_E/10**9,2),
                               np.round( self.aver_G/10**9,2),
                               np.round( self.aver_Nu,4),
                               np.round(self.aver_sum_cndVn_b6,8),
                               np.round(self.aver_Ty0/10**6,2),
                               np.round(self.aver_dEb*self.J2eV,4),
                               np.round(self.aver_delta_ss/10**6/self.model.Taylor_fac,2),
                               np.round(self.aver_delta_ss/10**6,2),
                               
                               np.round(self.std_V*10**30,2),
                               np.round(self.std_b*10**10,4),
                               np.round( self.std_E/10**9,2),
                               np.round( self.std_G/10**9,2),
                               np.round( self.std_Nu,4),
                               np.round(self.std_sum_cndVn_b6,8),
                               np.round(self.std_Ty0/10**6,2),
                               np.round(self.std_dEb*self.J2eV,4),
                               np.round(self.std_delta_ss/10**6/self.model.Taylor_fac,2),
                               np.round(self.std_delta_ss/10**6,2),
                               ]
                #print(calc_values)
                to_append = {**dict(zip(element_header,self.model.compositions.values)),**dict(zip(calc_header,calc_values))}
                #print(to_append)
                calc_data_all = calc_data_all.append(to_append,ignore_index=True)
        self.calc_data_all = calc_data_all
    def writedata(self,foutname):
        outputfile = foutname 
        self.calc_data_all.to_csv(outputfile,index=False,sep='\t',float_format='%.2f')
        print('-------------------------------Data Saved to {}-------------------------------'.format(outputfile))





class ss_edge_model_batch:
# calculate solid solution strengthening contribution for FCC/BCC CCAs
# Edge dislocation models
# FCC model: Varvenne-Leyson-Ghazisaeidi-Curtin 2016: http://dx.doi.org/10.1016/j.actamat.2016.09.046
# BCC model: Maresca-Curtin 2019: https://doi.org/10.1016/j.actamat.2019.10.015
    def __init__(self,
                 dislocation_properties,
                 exp_conditions,
                 comp_elements,comp_pst,
                 elements_data,
                 structure):
        
        # dislocation properties, alpha, f1, and f2
        self.alpha = float(dislocation_properties[0])
        self.f_tau = float(dislocation_properties[1])
        self.f_dEb = float(dislocation_properties[2])
        
        # experiment conditions, T, strain rate
        if len(exp_conditions[0]) > 1:
            self.T_range = np.array(exp_conditions[0])
        else:
            self.T = np.array(exp_conditions[0][0])
        self.ep = float(exp_conditions[1])
        self.ep0 = 10**4 #reference strain rate (/s)
        
        # some constants 
        self.boltzmann_J = 1.38064852*10**(-23) #J/K
        self.J2eV=6.2415093433*10**18 # covert J to eV
        
        # elemental data
        self.elements_order = comp_elements.columns.tolist()
        self.compositions = comp_elements #pandas df
        self.elements_data = copy.deepcopy(elements_data) #json
        
        # convert unit for properties
        # Vn: Å^3 to m^3
        # b: Å to m
        # a: Å to m
        # E: GPa to Pa
        # G: GPa to Pa
        
        for element_i in self.elements_order:
            self.elements_data[element_i]['Vn'] = elements_data[element_i]['Vn']*10**(-30)
            self.elements_data[element_i]['b'] = elements_data[element_i]['b']*10**(-10)
            self.elements_data[element_i]['a'] = elements_data[element_i]['a']*10**(-10)
            self.elements_data[element_i]['E'] = elements_data[element_i]['E']*10**(9)
            self.elements_data[element_i]['G'] = elements_data[element_i]['G']*10**(9)

        
        
        self.structure = structure

    def FCC_V_L_G_C_2016_analytical(self):
        # FCC model: Varvenne-Leyson-Ghazisaeidi-Curtin 2016: http://dx.doi.org/10.1016/j.actamat.2016.09.046
        
        self.prefac_ty0 = 0.051
        self.Taylor_fac = 3.06
        self.prefac_dEb = 0.274
        # averaged properties
        cn_Vn = []
        cn_nu = []
        cn_G = []
        cn_b = []
        cn_E = []
        for element_i in self.elements_order:
            cn_Vn.append(self.compositions[element_i]/100*self.elements_data[element_i]['Vn'])
            cn_nu.append(self.compositions[element_i]/100*self.elements_data[element_i]['nu'])
            cn_G.append(self.compositions[element_i]/100*self.elements_data[element_i]['G'])
            cn_b.append(self.compositions[element_i]/100*self.elements_data[element_i]['b'])
            cn_E.append(self.compositions[element_i]/100*self.elements_data[element_i]['E'])

        self.aver_E = sum(cn_E);
        self.aver_V = sum(cn_Vn);
        self.aver_G = sum(cn_G)
        self.aver_Nu = sum(cn_nu)
        self.aver_b = sum(cn_b)
        
        i = 0;cn_Delta_Vn2=[]
        for element_i in self.elements_order:
            cn_Delta_Vn2.append(self.compositions[element_i]/100*(self.elements_data[element_i]['Vn']-self.aver_V)**2)
            
        self.sum_cndVn_b6 = sum(cn_Delta_Vn2)/self.aver_b**6;
        q_nu = ((1 + self.aver_Nu)/(1 - self.aver_Nu))
        
        self.dEb = self.prefac_dEb * self.f_dEb * self.alpha**(1/3)  * self.aver_G * self.aver_b**3   * q_nu**(2/3) * self.sum_cndVn_b6**(1/3)
        self.Ty0 = self.prefac_ty0 * self.f_tau * self.alpha**(-1/3) * self.aver_G * q_nu**(4/3) * self.sum_cndVn_b6**(2/3)
        self.Ty0_pc = self.Taylor_fac * self.Ty0
        delta_ss_low_T = self.Ty0 * (1 - ((self.boltzmann_J*self.T)/(self.dEb) * np.log(self.ep0/self.ep))**(2/3) )
        delta_ss_high_T = self.Ty0 * np.exp(-1/0.57 * self.boltzmann_J*self.T/self.dEb * np.log(self.ep0/self.ep) )
        Ty_threshold = self.Ty0/2
        
        self.delta_ss = self.Taylor_fac*np.array([delta_ss_low_T[i] if delta_ss_low_T[i]>=Ty_threshold[i] else delta_ss_high_T[i] for i in range(len(Ty_threshold))])
        
        
    def BCC_M_C_2020_analytical(self):
        # BCC model: Maresca-Curtin-2019: https://doi.org/10.1016/j.actamat.2019.10.015
        
        self.prefac_ty0 = 0.051
        self.Taylor_fac = 3.06
        self.prefac_dEb = 0.274
        # averaged properties
        cn_Vn = []
        cn_nu = []
        cn_G = []
        cn_b = []
        cn_E = []
        for element_i in self.elements_order:
            cn_Vn.append(self.compositions[element_i]/100*self.elements_data[element_i]['Vn'])
            cn_nu.append(self.compositions[element_i]/100*self.elements_data[element_i]['nu'])
            cn_G.append(self.compositions[element_i]/100*self.elements_data[element_i]['G'])
            cn_b.append(self.compositions[element_i]/100*self.elements_data[element_i]['b'])
            cn_E.append(self.compositions[element_i]/100*self.elements_data[element_i]['E'])

        self.aver_E = sum(cn_E);
        self.aver_V = sum(cn_Vn);
        self.aver_G = sum(cn_G)
        self.aver_Nu = sum(cn_nu)
        self.aver_b = sum(cn_b)
        
        i = 0;cn_Delta_Vn2=[]
        for element_i in self.elements_order:
            cn_Delta_Vn2.append(self.compositions[element_i]/100*(self.elements_data[element_i]['Vn']-self.aver_V)**2)
            
        self.sum_cndVn_b6 = sum(cn_Delta_Vn2)/self.aver_b**6;
        q_nu = ((1 + self.aver_Nu)/(1 - self.aver_Nu))
        
        self.dEb = self.prefac_dEb * self.f_dEb * self.alpha**(1/3)  * self.aver_G * self.aver_b**3   * q_nu**(2/3) * self.sum_cndVn_b6**(1/3)
        self.Ty0 = self.prefac_ty0 * self.f_tau * self.alpha**(-1/3) * self.aver_G * q_nu**(4/3) * self.sum_cndVn_b6**(2/3)
        self.Ty0_pc = self.Taylor_fac * self.Ty0
        delta_ss_low_T = self.Ty0 * (1 - ((self.boltzmann_J*self.T)/(self.dEb) * np.log(self.ep0/self.ep))**(2/3) )
        delta_ss_high_T = self.Ty0 * np.exp(-1/0.57 * self.boltzmann_J*self.T/self.dEb * np.log(self.ep0/self.ep) )
        Ty_threshold = self.Ty0/2
        
        self.delta_ss = self.Taylor_fac*np.array([delta_ss_low_T[i] if delta_ss_low_T[i]>=Ty_threshold[i] else delta_ss_high_T[i] for i in range(len(Ty_threshold))])
        
    def calculate(self):
        self.calc_data_all = {}
        for T in self.T_range:
            self.T = T
            if self.structure == 'fcc':
                self.FCC_V_L_G_C_2016_analytical()
            elif self.structure == 'bcc':
                self.BCC_M_C_2020_analytical()
            self.calc_data = copy.deepcopy(comp_pst)
            self.calc_data['T'] = self.T*np.ones(len(self.calc_data))
            self.calc_data['V_ave'] = self.aver_V*10**30
            self.calc_data['b_ave'] = self.aver_b*10**10
            self.calc_data['E_ave'] = self.aver_E/10**9
            self.calc_data['G_ave'] = self.aver_G/10**9
            self.calc_data['nu_ave'] = self.aver_Nu
            self.calc_data['T'] = [self.T for i in range(len(self.Ty0))]
            self.calc_data['sum_cnVn^2_b6'] = [round(self.sum_cndVn_b6[i],8) for i in range(len(self.Ty0))]
            self.calc_data['Ty0'] = [round(self.Ty0[i]/10**6,2) for i in range(len(self.Ty0))]
            self.calc_data['Delta_Eb'] = [round(self.dEb[i]*self.J2eV,4) for i in range(len(self.Ty0))]
            self.calc_data['Delta_sigma_ss'] = [round(self.delta_ss[i]/10**6,2) for i in range(len(self.Ty0))]
            self.calc_data_all[str(self.T)] = self.calc_data
            
    def writedata(self,foutname):
        for key_i in self.calc_data_all.keys():
            calc_data = self.calc_data_all[key_i]
            outputfile = foutname+'_'+str(key_i)
            calc_data.to_csv(outputfile,index=False,sep='\t',float_format='%.2f')
            print('-------------------------------Data Saved to {}-------------------------------'.format(outputfile))






class ss_model_M_C_screw_pseudo_ternary:
    # BCC screw dislocation model: Maresca-Curtin 2019: https://doi.org/10.1016/j.actamat.2019.10.007
    # BCC_screw_Maresca-Curtin-2019
    # for pseudo-ternary prediction

    def __init__(self,
                 inputdata,
                 compositions,comp_pst
                ):

        # adjustable scalers
        self.kink_width = inputdata.adjustable_scalers['kink_width']  
        self.Delta_V_p_scaler = inputdata.adjustable_scalers['Delta_V_p_scaler']  
        self.Delta_E_p_scaler = inputdata.adjustable_scalers['Delta_E_p_scaler'] 
        self.comp_pst = comp_pst
        # some constants
        self.boltzmann_J = 1.38064852*10**(-23) #J/K
        self.boltzmann_eV = 8.617333262145e-5 #eV
        self.J2eV = 6.2415093433*10**18 # covert J to eV 
        self.eV2J = 1/self.J2eV
        
        # properties
        self.elements_order = compositions.columns.tolist()
        self.compositions = copy.deepcopy(compositions)
        self.elements_data = copy.deepcopy(inputdata.elements_data)
        cn_a = []
        cn_E_k = []
        cn_E_v = []
        cn_E_si = []
        cn_Delta_E_p = []
        cn_Delta_V_p = []
        
        for element_i in self.elements_order:
            cn_a.append(self.compositions[element_i]/100*self.elements_data[element_i]['a'])
            cn_E_k.append(self.compositions[element_i]/100*self.elements_data[element_i]['E_k'])
            cn_E_v.append(self.compositions[element_i]/100*self.elements_data[element_i]['E_f_v'])
            cn_E_si.append(self.compositions[element_i]/100*self.elements_data[element_i]['E_f_si'])
            cn_Delta_E_p.append(self.compositions[element_i]/100*self.elements_data[element_i]['Delta_E_p']**2)
            cn_Delta_V_p.append(self.compositions[element_i]/100*self.elements_data[element_i]['Delta_V_p'])

        self.a = sum(cn_a) * 10**(-10) 
        self.a_p = self.a*np.sqrt(2/3)  # Peierls spacing
        self.b = self.a*np.sqrt(3)/2    # burgers vector
        self.E_k = sum(cn_E_k) * self.eV2J 
        self.E_v = sum(cn_E_v) * self.eV2J 
        self.E_si = sum(cn_E_si) * self.eV2J 
        self.Delta_E_p = np.sqrt(sum(cn_Delta_E_p)) * self.Delta_E_p_scaler * self.eV2J 
        self.Delta_V_p = sum(cn_Delta_V_p) * self.Delta_E_p_scaler * self.eV2J /self.b
        
        # exp conditions
        self.T_range = np.arange(inputdata.conditions['temperature']['min'],
                           inputdata.conditions['temperature']['max']+inputdata.conditions['temperature']['inc'],
                           inputdata.conditions['temperature']['inc'])
        self.strain_r = inputdata.conditions['strain_r']  # strain rate
        self.strain_r_0 = 10**4                     # reference strain rate 10^4 /s
        
        self.w_k = self.kink_width * self.b                # kink width 
        self.xi_c = (1.083*self.E_k/self.Delta_E_p)**2*self.b                # characteristic length of dislocation segment 
        self.xi_si = self.xi_c * 15
        self.xi_v = self.xi_c * 7.5 
        
    def M_C_screw_model(self):
        
        # cross-kink
        # self-interstitial
        self.tau_xk_0_si = np.pi * self.E_si / (self.a_p * self.b * self.xi_si )
        self.tau_xk_si = self.tau_xk_0_si * (1-(self.Delta_H/self.E_si)**(2/3))
        # vacancy
        self.tau_xk_0_v = np.pi * self.E_v / (self.a_p * self.b * self.xi_v )
        self.tau_xk_v = self.tau_xk_0_v * (1-(self.Delta_H/self.E_v)**(2/3))
        # select the larger value from si or vacancy strengthening
        self.tau_xk_T = np.maximum(self.tau_xk_si,self.tau_xk_v)
        
        
        # kink glide
        self.tau_b = 1.08 * self.E_k / (self.a_p * self.b * self.xi_c)
        self.tau_k_0 = 6.3 * self.Delta_E_p / (self.a_p * self.b**2 * np.sqrt(self.w_k/self.b)) + self.tau_b
        self.Delta_E_k_0 = 1.37 * np.sqrt(self.w_k/self.b) * self.Delta_E_p 
        self.tau_k_low_T = self.tau_b + \
                           (self.tau_k_0 - self.tau_b) / \
                           (np.exp(0.89*self.Delta_H/self.Delta_E_k_0 + \
                                   0.5*(self.Delta_H/self.Delta_E_k_0)**(1/4) + 0.6)-1)
        
        self.tau_k_high_T = self.tau_b - \
                            (self.tau_k_0 - self.tau_b) * self.w_k / (5.75 * self.xi_c) * \
                            (self.Delta_H/self.Delta_E_k_0 - np.log(5.75*self.xi_c/self.w_k+1))
        
        self.tau_k_T = np.array([self.tau_k_low_T[i] if (self.tau_k_low_T[i]-self.tau_b[i])/(self.tau_k_0[i] - self.tau_b[i])>= 
                                 (1/(5.75 * self.xi_c[i]/self.w_k[i] + 1)) else 
                                 self.tau_k_high_T[i] for i in range(len(self.a))])
        
        # Peierls
        self.Delta_E_b_p =  (10*self.Delta_V_p*self.xi_c + 0.7 * self.E_k)**3/\
                            (20*self.Delta_V_p*self.xi_c + 0.7 * self.E_k)**2
        self.tau_p_0 = np.pi*self.Delta_V_p/(self.b*self.a_p ) + \
                       0.44 * self.E_k / (self.b*self.a_p * self.xi_c) * \
                      ( 1 - 5* self.Delta_V_p*self.xi_c/(20*self.Delta_V_p*self.xi_c+0.7*self.E_k))
        
        self.tau_p_T = self.tau_p_0 * (1-(self.Delta_H/self.Delta_E_b_p)**(2/3))
        
        # min of Peierls and kink glide
        self.min_tau_k_tau_p_T = np.minimum(self.tau_p_T,self.tau_k_T)
        
        # total strength
        self.tau_tot_T = np.maximum(self.min_tau_k_tau_p_T,np.zeros(len(self.a))) + self.tau_xk_T
    
    def calculate(self):
        self.calc_data_all = {}
        for T in self.T_range:
            self.T = T
            self.Delta_H =  self.boltzmann_J * self.T * np.log(self.strain_r_0/self.strain_r) #activation enthalpy
            self.M_C_screw_model()
            self.calc_data = copy.deepcopy(self.comp_pst)
            self.calc_data['a'] = np.round(self.a*10**10,4)
            self.calc_data['b'] = np.round(self.b*10**10,4)
            self.calc_data['a_p'] = np.round(self.a_p*10**10,4)
            self.calc_data['T'] = [self.T for i in range(len(self.calc_data))]
            self.calc_data['tau_y'] = np.round(self.tau_tot_T/10**6,2)
            self.calc_data['tau_k'] = np.round(self.tau_k_T/10**6,2)
            self.calc_data['tau_xk'] = np.round(self.tau_xk_T/10**6,2)
            self.calc_data['tau_p'] = np.round(self.tau_p_T/10**6,2)
            self.calc_data['E_k'] = np.round(self.E_k*self.J2eV,4)
            self.calc_data['E_f_v'] = np.round(self.E_v*self.J2eV,4)
            self.calc_data['E_f_si'] = np.round(self.E_si*self.J2eV,4)
            self.calc_data['Delta_E_p'] = np.round(self.Delta_E_p*(self.J2eV),4)
            self.calc_data_all[str(self.T)] = self.calc_data
    def writedata(self,foutname):
        for key_i in self.calc_data_all.keys():
            calc_data = self.calc_data_all[key_i]
            outputfile = foutname+'_'+str(key_i)
            calc_data.to_csv(outputfile,index=False,sep='\t',float_format='%.2f')
            print('-------------------------------Data Saved to {}-------------------------------'.format(outputfile))



class ss_model_M_C_screw:
    # BCC screw dislocation model: Maresca-Curtin 2019: https://doi.org/10.1016/j.actamat.2019.10.007


    def __init__(self,
                 inputdata
                ):

        # adjustable scalers
        self.kink_width = inputdata.adjustable_scalers['kink_width']  
        self.Delta_V_p_scaler = inputdata.adjustable_scalers['Delta_V_p_scaler']  
        self.Delta_E_p_scaler = inputdata.adjustable_scalers['Delta_E_p_scaler'] 
        
        # some constants
        self.boltzmann_J = 1.38064852*10**(-23) #J/K
        self.boltzmann_eV = 8.617333262145e-5 #eV
        self.J2eV = 6.2415093433*10**18 # covert J to eV 
        self.eV2J = 1/self.J2eV
        
        # properties
        self.a = inputdata.properties['a'] * 10**(-10)             #m    # lattice constant
        self.a_p = self.a*np.sqrt(2/3)  # Peierls spacing
        self.b = self.a*np.sqrt(3)/2
        
        self.E_k = inputdata.properties['E_k'] * self.eV2J             # J   # kink formation energy
        self.Delta_E_p = self.Delta_E_p_scaler * inputdata.properties['Delta_E_p'] * self.eV2J  # J   # screw-solute interaction
        self.Delta_V_p = self.Delta_V_p_scaler * inputdata.properties['Delta_V_p'] * self.eV2J /self.b# J/b   # Peierls barrier
        
        self.E_si = inputdata.properties['E_f_si']   * self.eV2J         #J   # formation energy of self-interstitial
        self.E_v = inputdata.properties['E_f_v']     * self.eV2J         #J   # formation energy of vacancy 
        
    
        # exp conditions
        self.T = np.arange(inputdata.conditions['temperature']['min'],
                           inputdata.conditions['temperature']['max']+inputdata.conditions['temperature']['inc'],
                           inputdata.conditions['temperature']['inc'])
        self.strain_r = inputdata.conditions['strain_r']  # strain rate
        self.strain_r_0 = 10**4                     # reference strain rate 10^4 /s
        

        self.Delta_H =  self.boltzmann_J * self.T * np.log(self.strain_r_0/self.strain_r)               #activation enthalpy
        self.w_k = self.kink_width * self.b                # kink width 
        self.xi_c = (1.083*self.E_k/self.Delta_E_p)**2*self.b                # characteristic length of dislocation segment 
        self.xi_si = self.xi_c * 15
        self.xi_v = self.xi_c * 7.5 
    def M_C_screw_model(self):
        
        # cross-kink
        # self-interstitial
        self.tau_xk_0_si = np.pi * self.E_si / (self.a_p * self.b * self.xi_si )
        self.tau_xk_si = self.tau_xk_0_si * (1-(self.Delta_H/self.E_si)**(2/3))
        # vacancy
        self.tau_xk_0_v = np.pi * self.E_v / (self.a_p * self.b * self.xi_v )
        self.tau_xk_v = self.tau_xk_0_v * (1-(self.Delta_H/self.E_v)**(2/3))
        # select the larger value from si or vacancy strengthening
        self.tau_xk_T = np.array([self.tau_xk_si[i] if self.tau_xk_si[i]>=self.tau_xk_v[i] else 
                                  self.tau_xk_v[i] for i in range(len(self.T)) ])
        
        
        # kink glide
        self.tau_b = 1.08 * self.E_k / (self.a_p * self.b * self.xi_c)
        self.tau_k_0 = 6.3 * self.Delta_E_p / (self.a_p * self.b**2 * np.sqrt(self.w_k/self.b)) + self.tau_b
        self.Delta_E_k_0 = 1.37 * np.sqrt(self.w_k/self.b) * self.Delta_E_p 
        self.tau_k_low_T = self.tau_b + \
                           (self.tau_k_0 - self.tau_b) / \
                           (np.exp(0.89*self.Delta_H/self.Delta_E_k_0 + \
                                   0.5*(self.Delta_H/self.Delta_E_k_0)**(1/4) + 0.6)-1)
        
        self.tau_k_high_T = self.tau_b - \
                            (self.tau_k_0 - self.tau_b) * self.w_k / (5.75 * self.xi_c) * \
                            (self.Delta_H/self.Delta_E_k_0 - np.log(5.75*self.xi_c/self.w_k+1))
        
        self.tau_k_T = np.array([self.tau_k_low_T[i] if (self.tau_k_low_T[i]-self.tau_b)/(self.tau_k_0 - self.tau_b)>= 
                                 (1/(5.75 * self.xi_c/self.w_k + 1)) else 
                                 self.tau_k_high_T[i] for i in range(len(self.T))])
        
        # Peierls
        self.Delta_E_b_p =  (10*self.Delta_V_p*self.xi_c + 0.7 * self.E_k)**3/\
                            (20*self.Delta_V_p*self.xi_c + 0.7 * self.E_k)**2
        self.tau_p_0 = np.pi*self.Delta_V_p/(self.b*self.a_p ) + \
                       0.44 * self.E_k / (self.b*self.a_p * self.xi_c) * \
                      ( 1 - 5* self.Delta_V_p*self.xi_c/(20*self.Delta_V_p*self.xi_c+0.7*self.E_k))
        
        self.tau_p_T = self.tau_p_0 * (1-(self.Delta_H/self.Delta_E_b_p)**(2/3))
        
        # min of Peierls and kink glide
        self.min_tau_k_tau_p_T = np.minimum(self.tau_p_T,self.tau_k_T)
        
        # total strength
        self.tau_tot_T = np.maximum(self.min_tau_k_tau_p_T,np.zeros(len(self.T))) + self.tau_xk_T

    def calculate(self):
        self.M_C_screw_model()
        self.calc_data = pd.DataFrame(data={})
        self.calc_data['T'] = self.T
        self.calc_data['tau_y'] = np.round(self.tau_tot_T/1e6,2)
    def writedata(self,foutname):
        self.calc_data.to_csv(foutname,index=False,sep='\t',float_format='%.2f')
        print('-------------------------------Data Saved to {}-------------------------------'.format(foutname))


        

class Suzuki_model_RWASM_jog_T:
    
    def __init__(self,
                adjustable_scalers,
                conditions,
                comp_elements,
                elements_data,
                T_l = None):
        
        # 
        self.elements_data = elements_data
        # conditions
        self.strain_r = float(conditions[1])
        self.T_range = np.array(conditions[0])
        self.elements_order = comp_elements.columns.tolist()
        self.compositions = comp_elements #pandas df
        # constants
        self.boltzmann_J = 1.380649e-23
        self.boltzmann_eV = 8.617333262145e-5
        self.J2eV = self.boltzmann_eV/self.boltzmann_J
        self.eV2J = 1/self.J2eV
        self.Debye = 5 * 10**(12) # Debye frequency /s
        self.T_l = T_l
        #adjustables
        self.rho = adjustable_scalers['dislocation_density']
        self.tau_i_exponent = adjustable_scalers['tau_i_exponent']
        self.trial_kappa_range = np.arange(adjustable_scalers['trial_kappa']['min'],
                                          adjustable_scalers['trial_kappa']['max']+adjustable_scalers['trial_kappa']['inc'],
                                          adjustable_scalers['trial_kappa']['inc'])
        self.trial_tau_k = adjustable_scalers['trial_tau_k'] * 1e6
        self.kink_width = adjustable_scalers['kink_width']

    def L(self,kappa_i):
        y = 1- stats.norm.cdf(kappa_i)
        return self.b/(3*y*self.c) 
    
    def tau_j(self,kappa_i):
        
        try: 
            liquidus_T = self.T_l[self.composition_idx]
        except:
            # if cannot read self.T_l[self.composition_idx] for that composition
            # assign a very high liquidus_T -> this will make jog mode impossible
            liquidus_T = 1e8
        Q = 19.12*self.boltzmann_J*liquidus_T
        D = 14*self.Debye*self.b**2 * np.exp(-Q/self.kT)
        time = 0.002/self.strain_r
        if np.sqrt(D*time) >= self.b:
            self.jog = True
            st = self.strain_r/(2.598*np.pi*self.rho*D)
            tau_j = (self.kT/(2*self.L(kappa_i)*self.b**2)) * np.log(st + np.sqrt(1+st**2) )
        else:
            self.jog = False
            tau_j = (self.E_int + self.E_vac)/(4*self.b*self.L(kappa_i))

        return tau_j
    
    def tau_y_k_j(self,tau_k_trial,kappa_i_list):
        
        self.Delta_V = lambda tau_k_trial,kappa_i: 3 * kappa_i**2 * self.E_w_i**2 * self.c / (2*tau_k_trial**2*self.a_p*self.b**2) + \
                                     tau_k_trial**2 * self.a_p**3 * self.b**4 * self.lambda_k**2 / (6*kappa_i**2 * self.E_w_i**2 * self.c)
        self.S = lambda tau_k_trial,kappa_i: 18 * kappa_i**2 * self.E_w_i**2 * self.c *self.kT /(self.a_p**3 * self.b**4 * self.lambda_k**2) * \
                 np.log( (5*np.pi*self.kT)**2 * self.Debye * self.a_p * self.b /((self.G*self.b*self.Delta_V(tau_k_trial,kappa_i))**2 * self.strain_r) )
        self.R = lambda kappa_i: 27 * kappa_i**4 * self.E_w_i**4 * self.c**2 / (self.a_p**4 * self.b**6 * self.lambda_k**2)
        # x[0] = tau_k_trial
        # x[1] = kappa_i
        self.tau_k_func = lambda tau_k_trial,kappa_i: tau_k_trial**4 + tau_k_trial*self.S(tau_k_trial,kappa_i) - self.R(kappa_i) 
        self.tau_k_kappa = []
        self.tau_j_kappa = []
        self.tau_y_kappa = []
        for kappa_i in kappa_i_list:
            self.res = optimize.root(self.tau_k_func, tau_k_trial, args=(kappa_i,))
            self.tau_k_kappa.append(self.res.x[0])
            self.tau_j_kappa.append(self.tau_j(kappa_i))
            self.tau_y_kappa.append(self.res.x[0] + self.tau_j(kappa_i))
        #print(self.tau_y_kappa)
        
    
    def phenomelogical_model_tau_y(self): 
        # tau_y = ( sum( tau_y_i**(1/q) ) )**q
        self.tau_y_tot = sum(self.tau_y_i**(1/self.tau_i_exponent))**self.tau_i_exponent
        
    def calculate(self):
        calc_data_all = []
        for composition_idx in range(len(self.compositions)):
            self.composition_idx = composition_idx
            tau_y_tot_T = []
            tau_y_i_T_list = []
            tau_k_i_T_list = []
            tau_j_i_T_list = []
            self.elements_kappa_i_convergence_record = pd.DataFrame(data={})
            for element_symbol in self.elements_data:
                self.elements_kappa_i_convergence_record[element_symbol] = {}
            jog_activated_T = []
            T_liquidus = []
            for T in self.T_range:
                self.T = T
                self.kT = self.boltzmann_J * self.T
                # record tau_y for every element
                tau_y_i = []
                tau_k_i = []
                tau_j_i = []

                for element_symbol in self.elements_data:
                    element_i = self.elements_data[element_symbol]
                    #print(element_i)
                    # calculate the yield strength contribution for every element
                    # according to concentration
                    # setup properties for every element
                    self.E_f_v = element_i['E_f_v'] * self.eV2J #J
                    self.E_f_si = element_i['E_f_si'] * self.eV2J # J
                    self.a_0 = element_i['a']*1e-10#element_i['a_0'] * 10**(-10) # unit: m
                    self.E_w_i = element_i['E_w_i'] * self.eV2J#element_i['E_w'] * self.eV2J # J
                    self.c = self.compositions[element_symbol][composition_idx]/100
                    self.G = element_i['G'] * 10**9 # Pa
                    self.nu = element_i['nu']
                    self.b = self.a_0 * np.sqrt(3) / 2
                    self.a_p = self.a_0 * np.sqrt(2/3)
                    #self.E_vac = 0.6 * self.eV2J / 10**(-10) # test NbTiZr
                    #self.E_int = 0.9 * self.eV2J / 10**(-10) # test NbTiZr
                    self.E_vac = 0.707 * self.E_f_v  /self.b + self.G * self.b**2 / (np.pi*(1-self.nu)) * np.log(1.5)
                    self.E_int = 0.707 * self.E_f_si /self.b + self.G * self.b**2 / (np.pi*(1-self.nu)) * np.log(1.5)
                    self.lambda_k = self.b * self.kink_width

                    # record the optimization results for post-processing
                    tau_k_list = []
                    tau_j_list = []
                    tau_y_list = []
                    optimized_kappa_list = []

                    # start to optimize tau_k for every trial kappa

                    self.tau_y_k_j(self.trial_tau_k, self.trial_kappa_range)
                    tau_k_list = np.array(self.tau_k_kappa)/1e6
                    tau_j_list = np.array(self.tau_j_kappa)/1e6
                    tau_y_list = np.array(self.tau_y_kappa)/1e6
                    kappa_list = self.trial_kappa_range

                    # optimize tau_y over kappa, this finds the true tau_y for each element
                    kappa_sort, tau_y_sort, tau_j_sort, tau_k_sort  = zip(*sorted(zip(kappa_list, tau_y_list, tau_j_list, tau_k_list)))
                    ###
                    polyfit = np.polyfit(kappa_sort,tau_y_sort,4)
                    npfit = np.poly1d(polyfit)
                    guess_kappa = (self.trial_kappa_range[0]+self.trial_kappa_range[1])/2
                    res = optimize.fmin_slsqp(npfit,guess_kappa,
                                              bounds=([(self.trial_kappa_range[0],self.trial_kappa_range[-1])]),
                                              iprint=-1)
                    optimized_kappa = res[0]
                    tau_y_i.append(npfit(optimized_kappa))
                    tau_k_i.append(npfit(optimized_kappa)-(self.tau_j(optimized_kappa)/1e6))
                    tau_j_i.append(self.tau_j(optimized_kappa)/1e6)

                    # record for convergence check
                    self.elements_kappa_i_convergence_record[element_symbol]['kappa_'+str(self.T)] = None # strange thing here, only by setting None it records the first row of data
                    self.elements_kappa_i_convergence_record[element_symbol]['tau_y_'+str(self.T)] = None
                    self.elements_kappa_i_convergence_record[element_symbol]['kappa_'+str(self.T)] = kappa_sort
                    self.elements_kappa_i_convergence_record[element_symbol]['tau_y_'+str(self.T)] = tau_y_sort

                tau_y_i_T_list.append(tau_y_i)
                tau_k_i_T_list.append(tau_k_i)
                tau_j_i_T_list.append(tau_j_i)
                jog_activated_T.append(self.jog)
                try:
                    T_liquidus.append(self.T_l[composition_idx])
                except:
                    T_liquidus.append(1e8)
                self.tau_y_i = np.array(tau_y_i)
                self.phenomelogical_model_tau_y()
                tau_y_tot_T.append(self.tau_y_tot)
            self.tau_y_tot_T = np.array(tau_y_tot_T)
            self.tau_y_i_T_list = np.array(tau_y_i_T_list).transpose()
            self.tau_k_i_T_list = np.array(tau_k_i_T_list).transpose()
            self.tau_j_i_T_list = np.array(tau_j_i_T_list).transpose()
            self.jog_activated_T = (jog_activated_T)
            self.T_liquidus = np.array(T_liquidus)
            self.T_liquidus[self.T_liquidus==1e8]=np.nan
            # record data
            calc_data_comp_idx = {}
            for element_symbol in self.elements_data:
                calc_data_comp_idx[element_symbol] = np.ones(len(self.T_range))*self.compositions[element_symbol][composition_idx]
            calc_data_comp_idx["T"] = self.T_range
            calc_data_comp_idx["tau_y"] = np.round(self.tau_y_tot_T,2)
            for i, element_symbol in zip(range(len(self.elements_data)),self.elements_data):
                calc_data_comp_idx["tau_y_"+str(element_symbol)] = np.round(self.tau_y_i_T_list[i],2)
                calc_data_comp_idx["tau_k_"+str(element_symbol)] = np.round(self.tau_k_i_T_list[i],2)
                calc_data_comp_idx["tau_j_"+str(element_symbol)] = np.round(self.tau_j_i_T_list[i],2)
            calc_data_comp_idx['jog_activated'] = self.jog_activated_T
            calc_data_comp_idx['T_liquidus'] = self.T_liquidus
            calc_data_all.append(pd.DataFrame(data=calc_data_comp_idx))
        self.calc_data_all = pd.concat(calc_data_all)
    def writedata(self,foutname):
        outputfile = foutname 
        self.calc_data_all.to_csv(outputfile,index=False,sep='\t',float_format='%.2f')
        print('-------------------------------Data Saved to {}-------------------------------'.format(outputfile))






class Suzuki_model_RWASM_ternary:
    
    def __init__(self,
                inputdata,
                comp_elements,comp_pst,
                T_l = None):
        
        
        self.element_composition = comp_elements
        self.elements_data = inputdata.elements_data
        self.comp_pst = comp_pst
        if T_l is None:
            print('no Tl (liquidus temperature) supplied, cannot predict high-temperature jog-dragging mode.')
            self.T_l = np.ones(len(self.element_composition))*1e8
        else:
            self.T_l = T_l
        # 
        
        # conditions
        self.strain_r = inputdata.conditions['strain_r']
        self.T_range = np.arange(inputdata.conditions['temperature']['min'],
                           inputdata.conditions['temperature']['max']+inputdata.conditions['temperature']['inc'],
                           inputdata.conditions['temperature']['inc'])
        
        # constants
        self.boltzmann_J = 1.380649e-23
        self.boltzmann_eV = 8.617333262145e-5
        self.J2eV = self.boltzmann_eV/self.boltzmann_J
        self.eV2J = 1/self.J2eV
        self.Debye = 5 * 10**(12) # Debye frequency /s

        #adjustables
        self.rho = inputdata.adjustable_scalers['dislocation_density']
        self.tau_i_exponent = inputdata.adjustable_scalers['tau_i_exponent']
        self.trial_kappa_range = np.arange(inputdata.adjustable_scalers['trial_kappa']['min'],
                                          inputdata.adjustable_scalers['trial_kappa']['max']+inputdata.adjustable_scalers['trial_kappa']['inc'],
                                          inputdata.adjustable_scalers['trial_kappa']['inc'])
        self.trial_tau_k = inputdata.adjustable_scalers['trial_tau_k'] * 1e6
        self.kink_width = inputdata.adjustable_scalers['kink_width']

        
        
    def L(self,kappa_i):
        f = lambda x: np.exp(-x**2/2)/np.sqrt(2*np.pi)
        y = integrate.quad(f,kappa_i,np.inf)
        return self.b/(3*y[0]*self.c) 
        
    def tau_j(self,kappa_i):
        
        try: 
            liquidus_T = self.T_l[self.composition_idx]
        except:
            # if cannot read self.T_l[self.composition_idx] for that composition
            # assign a very high liquidus_T -> this will make jog mode impossible
            liquidus_T = 1e8
        Q = 19.12*self.boltzmann_J*liquidus_T
        D = 14*self.Debye*self.b**2 * np.exp(-Q/self.kT)
        time = 0.002/self.strain_r
        if np.sqrt(D*time) >= self.b:
            self.jog = True
            st = self.strain_r/(2.598*np.pi*self.rho*D)
            tau_j = (self.kT/(2*self.L(kappa_i)*self.b**2)) * np.log(st + np.sqrt(1+st**2) )
        else:
            self.jog = False
            tau_j = (self.E_int + self.E_vac)/(4*self.b*self.L(kappa_i))

        return tau_j
    
    def tau_y_k_j(self,tau_k_trial,kappa_i_list):

        self.Delta_V = lambda tau_k_trial,kappa_i: 3 * kappa_i**2 * self.E_w**2 * self.c / (2*tau_k_trial**2*self.a_p*self.b**2) + \
                                     tau_k_trial**2 * self.a_p**3 * self.b**4 * self.lambda_k**2 / (6*kappa_i**2 * self.E_w**2 * self.c)
        self.S = lambda tau_k_trial,kappa_i: 18 * kappa_i**2 * self.E_w**2 * self.c *self.kT /(self.a_p**3 * self.b**4 * self.lambda_k**2) * \
                 np.log( (5*np.pi*self.kT)**2 * self.Debye * self.a_p * self.b /((self.G*self.b*self.Delta_V(tau_k_trial,kappa_i))**2 * self.strain_r) )
        self.R = lambda kappa_i: 27 * kappa_i**4 * self.E_w**4 * self.c**2 / (self.a_p**4 * self.b**6 * self.lambda_k**2)

        self.tau_k_func = lambda tau_k_trial,kappa_i: tau_k_trial**4 + tau_k_trial*self.S(tau_k_trial,kappa_i) - self.R(kappa_i) 
        self.tau_k_kappa = []
        self.tau_j_kappa = []
        self.tau_y_kappa = []
        for kappa_i in kappa_i_list:
            self.res = optimize.root(self.tau_k_func, tau_k_trial, args=(kappa_i,))
            self.tau_k_kappa.append(self.res.x[0])
            self.tau_j_kappa.append(self.tau_j(kappa_i))
            self.tau_y_kappa.append(self.res.x[0] + self.tau_j(kappa_i))
        
    
    def phenomelogical_model_tau_y(self): 
        # tau_y = ( sum( tau_y_i**(1/q) ) )**q
        self.tau_y_tot = [sum(tau_y_i**(1/self.tau_i_exponent))**self.tau_i_exponent for tau_y_i in self.tau_y_i_pst]
        
    def calculate(self):
        self.calc_data_all = {}
        for T in self.T_range:
            self.T = T
            self.kT = self.boltzmann_J * self.T
            self.tau_y_i_pst = [] # record all compositions shape(len(composition),len(element))
            tau_y_i_pst = []
            self.jog_i = []
            for i in range(len(self.element_composition)):
                self.composition_idx = i
                # record tau_y for every element at composition i
                tau_y_i = []
                for element_symbol in self.element_composition.columns:
                    element_i = self.elements_data[element_symbol]
                    #print(element_i)
                    # calculate the yield strength contribution for every element
                    # according to concentration
                    # setup properties for every element
                    self.c = self.element_composition[element_symbol][i]/100
                    if self.c == 0:
                        tau_y_i.append(0)
                        continue
                    self.E_f_v = element_i['E_f_v'] * self.eV2J #J
                    self.E_f_si = element_i['E_f_si'] * self.eV2J # J
                    self.a_0 = element_i['a']*1e-10#element_i['a_0'] * 10**(-10) # unit: m
                    self.E_w = element_i['E_w_i'] * self.eV2J#element_i['E_w_i'] * self.eV2J # J

                    self.G = element_i['G'] * 10**9 # Pa
                    self.nu = element_i['nu']
                    self.b = self.a_0 * np.sqrt(3) / 2
                    self.a_p = self.a_0 * np.sqrt(2/3)
                    #self.E_vac = 0.6 * self.eV2J / 10**(-10) # test NbTiZr
                    #self.E_int = 0.9 * self.eV2J / 10**(-10) # test NbTiZr
                    self.E_vac = 0.707 * self.E_f_v  /self.b + self.G * self.b**2 / (np.pi*(1-self.nu)) * np.log(1.5)
                    self.E_int = 0.707 * self.E_f_si /self.b + self.G * self.b**2 / (np.pi*(1-self.nu)) * np.log(1.5)
                    self.lambda_k = self.b * self.kink_width

                    # record the optimization results for post-processing
                    tau_k_list = []
                    tau_j_list = []
                    tau_y_list = []
                    optimized_kappa_list = []

                    # start to optimize tau_k for every trial kappa
                    self.tau_y_k_j(self.trial_tau_k, self.trial_kappa_range)
                    tau_k_list = np.array(self.tau_k_kappa)/1e6
                    tau_j_list = np.array(self.tau_j_kappa)/1e6
                    tau_y_list = np.array(self.tau_y_kappa)/1e6
                    kappa_list = self.trial_kappa_range

                    # optimize tau_y over kappa, this finds the true tau_y for each element
                    kappa_sort, tau_y_sort, tau_j_sort, tau_k_sort  = zip(*sorted(zip(kappa_list, tau_y_list, tau_j_list, tau_k_list)))
                    ###
                    polyfit = np.polyfit(kappa_sort,tau_y_sort,4)
                    npfit = np.poly1d(polyfit)
                    guess_kappa = (self.trial_kappa_range[0]+self.trial_kappa_range[1])/2
                    res = optimize.fmin_slsqp(npfit,guess_kappa,
                                              bounds=([(self.trial_kappa_range[0],self.trial_kappa_range[-1])]),
                                              iprint=-1)
                    optimized_kappa = res[0]
                    tau_y_i.append(npfit(optimized_kappa))

                self.tau_y_i = np.array(tau_y_i)

                self.jog_i.append(self.jog)
                
                tau_y_i_pst.append(self.tau_y_i)
            self.tau_y_i_pst = np.array(tau_y_i_pst)
            self.phenomelogical_model_tau_y()
            # record data for temperature
            self.calc_data = copy.deepcopy(self.comp_pst)
            for idx in range(len(self.element_composition.columns)):
                # elemental contribution to tau_y
                self.calc_data['tau_y_'+self.element_composition.columns[idx]] = np.round(self.tau_y_i_pst.transpose()[idx],2)
            self.calc_data['T'] = np.ones(len(self.calc_data)) * self.T
            self.calc_data['tau_y'] = np.round(self.tau_y_tot,2)
            self.calc_data['T_liquidus'] = self.T_l
            self.calc_data[self.calc_data['T_liquidus']==1e8]['T_liquidus']=np.nan
            self.calc_data['Jog_mode'] = self.jog_i
            self.calc_data_all[str(self.T)] = self.calc_data

    def writedata(self,foutname):
        for key_i in self.calc_data_all.keys():
            calc_data = self.calc_data_all[key_i]
            outputfile = foutname+'_'+str(key_i)
            calc_data.to_csv(outputfile,index=False,sep='\t',float_format='%.2f')
            print('-------------------------------Data Saved to {}-------------------------------'.format(outputfile))



class ss_model_M_C_screw_single:
    # BCC screw dislocation model: Maresca-Curtin 2019: https://doi.org/10.1016/j.actamat.2019.10.007
    # similar to pseudo-ternary, but for single-composition calculations using c-weighted properties.

    def __init__(self,
                 adjustable_scalers,
                  conditions, # conditions = [[list of T],strain_r]
                  compositions,  # dict
                  elements_data   # dict 
                ):
        
        
        # adjustable scalers
        self.kink_width = adjustable_scalers['kink_width']  
        self.Delta_V_p_scaler = adjustable_scalers['Delta_V_p_scaler']  
        self.Delta_E_p_scaler = adjustable_scalers['Delta_E_p_scaler'] 
        
        # some constants
        self.boltzmann_J = 1.38064852*10**(-23) #J/K
        self.boltzmann_eV = 8.617333262145e-5 #eV
        self.J2eV = self.boltzmann_eV/self.boltzmann_J # covert J to eV 
        self.eV2J = 1/self.J2eV
        

        # properties
        self.elements_order = compositions.columns.to_list()
        self.compositions = copy.deepcopy(compositions)
        self.elements_data = copy.deepcopy(elements_data)
        cn_a = []
        cn_E_k = []
        cn_E_v = []
        cn_E_si = []
        cn_Delta_E_p = []
        cn_Delta_V_p = []
        comp_normalizer = sum([self.compositions[element_i].values for element_i in self.elements_order])
        for element_i in self.elements_order:
            cn_a.append(self.compositions[element_i].values/comp_normalizer*self.elements_data[element_i]['a'])
            cn_E_k.append(self.compositions[element_i].values/comp_normalizer*self.elements_data[element_i]['E_k'])
            cn_E_v.append(self.compositions[element_i].values/comp_normalizer*self.elements_data[element_i]['E_f_v'])
            cn_E_si.append(self.compositions[element_i].values/comp_normalizer*self.elements_data[element_i]['E_f_si'])
            cn_Delta_E_p.append(self.compositions[element_i].values/comp_normalizer*self.elements_data[element_i]['Delta_E_p']**2)
            cn_Delta_V_p.append(self.compositions[element_i].values/comp_normalizer*self.elements_data[element_i]['Delta_V_p'])

        self.a = sum(cn_a) * 10**(-10) 
        self.a_p = self.a*np.sqrt(2/3)  # Peierls spacing
        self.b = self.a*np.sqrt(3)/2    # burgers vector
        self.E_k = sum(cn_E_k) * self.eV2J 
        self.E_v = sum(cn_E_v) * self.eV2J 
        self.E_si = sum(cn_E_si) * self.eV2J 
        self.Delta_E_p = np.sqrt(sum(cn_Delta_E_p)) * self.Delta_E_p_scaler * self.eV2J 
        self.Delta_V_p = sum(cn_Delta_V_p) * self.Delta_V_p_scaler * self.eV2J /self.b
        
        # exp conditions
        self.T_range = conditions[0]
        self.strain_r =   conditions[1]  # strain rate
        self.strain_r_0 = 1e4                     # reference strain rate 10^4 /s
        
        self.w_k = self.kink_width * self.b                # kink width 
        self.xi_c = (1.083*self.E_k/self.Delta_E_p)**2*self.b     # characteristic length of dislocation segment 
        self.xi_si = self.xi_c * 15
        self.xi_v = self.xi_c * 7.5 
        
    def M_C_screw_model(self):
        
        # cross-kink
        # self-interstitial
        self.Delta_H =  self.boltzmann_J * self.T * np.log(self.strain_r_0/self.strain_r)               #activation enthalpy
        self.tau_xk_0_si = np.pi * self.E_si / (self.a_p * self.b * self.xi_si )
        self.tau_xk_si = self.tau_xk_0_si * (1-(self.Delta_H/self.E_si)**(2/3))
        # vacancy
        self.tau_xk_0_v = np.pi * self.E_v / (self.a_p * self.b * self.xi_v )
        self.tau_xk_v = self.tau_xk_0_v * (1-(self.Delta_H/self.E_v)**(2/3))
        
        # select the larger value from si or vacancy strengthening
        
        self.tau_xk_T = np.maximum(self.tau_xk_si,self.tau_xk_v)
        
        
                                
        
        
        # kink glide
        self.tau_b = 1.08 * self.E_k / (self.a_p * self.b * self.xi_c)
        self.tau_k_0 = 6.3 * self.Delta_E_p / (self.a_p * self.b**2 * np.sqrt(self.w_k/self.b)) + self.tau_b
        self.Delta_E_k_0 = 1.37 * np.sqrt(self.w_k/self.b) * self.Delta_E_p 
        self.tau_k_low_T = self.tau_b + \
                           (self.tau_k_0 - self.tau_b) / \
                           (np.exp(0.89*self.Delta_H/self.Delta_E_k_0 + \
                                   0.5*(self.Delta_H/self.Delta_E_k_0)**(1/4) + 0.6)-1)
        
        self.tau_k_high_T = self.tau_b - \
                            (self.tau_k_0 - self.tau_b) * self.w_k / (5.75 * self.xi_c) * \
                            (self.Delta_H/self.Delta_E_k_0 - np.log(5.75*self.xi_c/self.w_k+1))
        if (self.tau_k_low_T-self.tau_b)/(self.tau_k_0 - self.tau_b)>=(1/(5.75 * self.xi_c/self.w_k + 1)):
            self.tau_k_T = self.tau_k_low_T
        else:
            self.tau_k_T = self.tau_k_high_T
        
        # Peierls
        self.Delta_E_b_p =  (10*self.Delta_V_p*self.xi_c + 0.7 * self.E_k)**3/\
                            (20*self.Delta_V_p*self.xi_c + 0.7 * self.E_k)**2
        self.tau_p_0 = np.pi*self.Delta_V_p/(self.b*self.a_p ) + \
                       0.44 * self.E_k / (self.b*self.a_p * self.xi_c) * \
                      ( 1 - 5* self.Delta_V_p*self.xi_c/(20*self.Delta_V_p*self.xi_c+0.7*self.E_k))
        
        self.tau_p_T = self.tau_p_0 * (1-(self.Delta_H/self.Delta_E_b_p)**(2/3))
        
        # min of Peierls and kink glide
        self.min_tau_k_tau_p_T = np.minimum(self.tau_p_T,self.tau_k_T)
        
        # total strength
        self.tau_tot_T = np.maximum(self.min_tau_k_tau_p_T,0) + self.tau_xk_T
        self.V_act = self.a_p * self.b * self.w_k / 2.3 * ((self.tau_k_0-self.tau_b)/(self.tau_tot_T-self.tau_b)-1)

    def calculate(self):
        calc_data_all = []
        for composition_idx in range(len(self.compositions)):
            self.composition_idx = composition_idx
            for T in self.T_range:
                calc_data = self.compositions.iloc[self.composition_idx].to_dict()
                self.T = T
                self.M_C_screw_model()
                calc_data['T'] = self.T
                calc_data['tau_y'] = np.round(self.tau_tot_T/1e6,2)
                calc_data['tau_k'] = np.round(self.tau_k_T/1e6,2)
                calc_data['tau_xk'] = np.round(self.tau_xk_T/1e6,2)
                calc_data['tau_p'] = np.round(self.tau_p_T/1e6,2)
                calc_data['tau_xk_v'] = np.round(self.tau_xk_v/1e6,2)
                calc_data['tau_xk_si'] = np.round(self.tau_xk_si/1e6,2)
                calc_data['V_act'] = np.round(self.V_act/self.b**3,2) 
                calc_data_all.append(pd.DataFrame(calc_data))
                
                
        self.calc_data_all = (pd.concat(calc_data_all,ignore_index=True))
    def writedata(self,foutname):
        self.calc_data_all.to_csv(foutname,index=False,sep='\t',float_format='%.2f')
        print('-------------------------------Data Saved to {}-------------------------------'.format(foutname))




