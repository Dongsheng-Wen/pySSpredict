import os
os.environ['TC20A_HOME'] = '/apps/cent7/thermocalc/2020a'
os.environ['LSHOST']='mooring.ecn.purdue.edu'
from tc_python import *
import numpy as np
import pandas as pd 

start_api_server()

class liquidus_wrapper:
    def __init__(self,comp_elements,tcdatabase="TCHEA4"):
        start_api_server()
        self.elements = [i for i in comp_elements.keys()]
        #print(self.elements)
        try:
            self.comp_elements = pd.DataFrame(data=comp_elements) # create pd.df from dict
        except:
            self.comp_elements = pd.DataFrame(data=comp_elements,index=[0]) 
        self.tcdatabase = tcdatabase
        self.system = SetUp().select_database_and_elements(tcdatabase,self.elements).get_system()
        
    def calc_liquidus(self):
        T_liquidus = []
        for i in range(len(self.comp_elements)):

            self.TC_calc = self.system.with_single_equilibrium_calculation()

            for element in self.elements[1:]: # must exclude the first elements, otherwise TC raises degree of freedom error
                print(element)
                self.TC_calc.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(element),self.comp_elements[element][i])

            self.TC_calc.remove_condition(ThermodynamicQuantity.temperature()) \
                        .set_phase_to_fixed('LIQUID', 1)
            TL = self.TC_calc.calculate().get_value_of(ThermodynamicQuantity.temperature())    
            T_liquidus.append(round(TL,2))

        self.T_liquidus = np.array(T_liquidus)


    def save_data(self):
        self.data = self.comp_elements
        self.data['T_liquidus'] = self.T_liquidus 
        fout_name = "".join([str(elem) for elem in self.elements])+'_{}_liquidus_T.csv'.format(self.tcdatabase)
        self.data.to_csv(fout_name,sep=',')

    def stop_server(self):
        stop_api_server() 


def main():
    df = data={
        "Nb": np.array([100,99,34]),
        "Ti": np.array([0,1,33]),
        "Zr": np.array([0,0,33])
        }
    lw = liquidus_wrapper(df)
    lw.calc_liquidus()
    lw.stop_server()
    print(lw.T_liquidus)



if __name__ == '__main__':
    main()

