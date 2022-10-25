import os
os.environ['TC21A_HOME'] = '/apps/cent7/thermocalc/2021a'
os.environ['LSHOST']='mooring.ecn.purdue.edu'
try:
    from tc_python import *
except:
    print('cannot import TCPython')
import numpy as np
import pandas as pd 
import json
import time 
try:
    start_api_server()
except:
    pass 
class ternary_isothermal_pd:
    # inspired by pyTCPlotter 
    # now should work for ternary alloys
    # need to expand to pseudo-ternary
    def __init__(self,elements,temperature,tcdatabase="TCHEA4"):
        start_api_server()
        self.temp = temperature 
        self.elements = elements
        # elements should be in the correct order for plotting the diagram 
        # ['A','B','C'] -> ['bottom right', 'top','bottom left']
        # this order is reflected in the output file name 
        #print(self.elements)
        self.tcdatabase = tcdatabase
        self.system = SetUp().select_database_and_elements(tcdatabase,self.elements).get_system()
        
    def calc_diagram(self):
        # most of the codes from pyTCPlotter make_ternary_diagram
        print('Set up phase diagram calculation at {}K'.format(self.temp))
        
        # setup the axes
        # self.elements = ['X', 'Y', 'Z']
        # order and positions in the phase diagram 
        # X - bottom right
        # Y - top
        # Z - bottom left 
        PD_calculator = self.system.with_phase_diagram_calculation()
        PD_calculator.with_first_axis(CalculationAxis(ThermodynamicQuantity.mole_fraction_of_a_component(self.elements[0])).
                                      set_min(0).
                                      set_max(1))
        PD_calculator.with_second_axis(CalculationAxis(ThermodynamicQuantity.mole_fraction_of_a_component(self.elements[1])).
                                       set_min(0).
                                       set_max(1))
        PD_calculator.set_condition(ThermodynamicQuantity.temperature(),self.temp)
        self.PD_result = PD_calculator.calculate()

    def make_diagram_data(self):
        # most of the codes from pyTCPlotter make_ternary_diagram

        PD_grouping = self.PD_result.get_values_grouped_by_stable_phases_of(
                                    ThermodynamicQuantity.mole_fraction_of_a_component(self.elements[0]),
                                    ThermodynamicQuantity.mole_fraction_of_a_component(self.elements[1]))
        # get the phase boundary lines
        self.pd_data = {'phase_boundaries':{} }
        for group in PD_grouping.get_lines().values():
            z = np.array([1-i-j for i,j in zip(group.x,group.y)])
            points=list(zip(100*np.array(group.x),100*np.array(group.y),100*z))
            phase_label = group.label
            self.pd_data['phase_boundaries'][phase_label] = points 

        # get the tielines
        z = np.array([1-i-j for i,j in zip(PD_grouping.get_tie_lines().x,PD_grouping.get_tie_lines().y)])
        x = np.array(PD_grouping.get_tie_lines().x)
        y = np.array(PD_grouping.get_tie_lines().y)
        tiel_points = list(zip(x*100,y*100,z*100))
        self.pd_data['tielines'] = tiel_points


    def save_data(self):

        fout_name = "".join([str(elem) for elem in self.elements])+'_{}_phase_diagram_{}K.json'.format(self.tcdatabase,self.temp)

        with open(fout_name, "w") as outfile:
            json.dump(self.pd_data, outfile)
        print('Saved phase diagram data to {}'.format(fout_name))
    def stop_server(self):
        stop_api_server() 


class single_point_isothermo: 
    # some codes copy from tc_python_plotter.py 

    def __init__(self,
        temperature = None,
        tcdatabase="TCHEA4"):

        start_api_server()
        self.temp = temperature
        self.tcdatabase = tcdatabase
    def set_system(self,Setup_TC,comp):
        # comp: dict; example: {"A":0.4, "B":0.4, "C":0.2}
        self.comp = comp 
        self.elements = [e for e in comp.keys()]
        
        self.system = Setup_TC.select_database_and_elements(self.tcdatabase,self.elements).get_system()
        self.calculator = self.system.with_single_equilibrium_calculation()
        self.calculator_1 = self.system.with_single_equilibrium_calculation()
        self.calculator_2 = self.system.with_single_equilibrium_calculation()
        self.calculator_3 = self.system.with_single_equilibrium_calculation()
        self.calculator_4 = self.system.with_single_equilibrium_calculation()
        self.calculator_5 = self.system.with_single_equilibrium_calculation()
        self.calculator_6 = self.system.with_single_equilibrium_calculation()
        self.calculator_7 = self.system.with_single_equilibrium_calculation()
        self.calculator_8 = self.system.with_single_equilibrium_calculation()
        self.calculator_9 = self.system.with_single_equilibrium_calculation()
                
    def set_tcdatabase(self,tcdatabase="TCHEA4"):
        self.tcdatabase = tcdatabase

    def set_comp(self,comp=None):
        # this is different from set_system()
        # this set composition only 
        # save lots of time loading the TC database
        
        if comp is not None: 
            self.comp = comp 
        i=1
        for ele in self.comp:
            if i < len(self.comp):
                self.calculator.set_condition(ThermodynamicQuantity.mole_fraction_of_a_component(ele), 
                    self.comp[ele])
                i+=1
            else:
                pass 
    def set_temperature(self,temperature):
        print('Set up single point calculation at {}K'.format(self.temp))
        self.temp = temperature # Kelvin 

    def calculate(self):
        # most of the codes from pyTCPlotter make_ternary_diagram
        self.stable_phases = {}
        try:
            self.calculator.set_condition(ThermodynamicQuantity.temperature(),self.temp)
            self.SE_result = self.calculator.calculate()
            phases = self.SE_result.get_stable_phases()
            for i in range(len(phases)):
                self.stable_phases[str(int(i+1))] = phases[i]
            #self.SE_result.invalidate()
        except:
            print('cannot calculate, pass...')
        self.calculator = self.system.with_single_equilibrium_calculation()

    def stop_server(self):
        stop_api_server() 