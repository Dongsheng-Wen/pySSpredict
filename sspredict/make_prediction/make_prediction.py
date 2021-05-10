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

# -----Dongsheng Wen-----
# -----2021-04-30--------
class make_prediction:

    def __init__(self,inputfilename):
        self.inputfilename = inputfilename
        self.m = master_read_input(self.inputfilename)
        self.m.assign_model_mode()

    def predict(self):
        # assign model to calculate based on tags
        if self.m.mode == "edge_ternary":
            print('---------Start Calculation for Ternary----------')
            print('-----------Use Edge Dislocation Model-----------')
            input_data = read_inputjson_edge_pseudo_ternary(self.inputfilename)
            input_data.check_integrity_curtin_edge()
            input_data.grab_properties_curtin_edge()
            self.input_data = input_data
            mesh = build_mesh_ternary(input_data.increment,
                          input_data.psA_range,input_data.psB_range,input_data.psC_range,
                          input_data.psA_elements,input_data.psB_elements,input_data.psC_elements,
                          input_data.psA_ratio,input_data.psB_ratio,input_data.psC_ratio)
            comp_tmp = mesh.make_mesh()
            comp_pst, comp_elements = mesh.comp_assign()
            if input_data.uncertainty_levels == [0,0]:
                self.model = ss_edge_model(input_data.dislocation_properties,input_data.exp_conditions,
                            comp_elements,comp_pst,input_data.element_data,input_data.structure)
            else: 
                self.model = ss_edge_model_w_uncertainty(ss_edge_model,input_data.dislocation_properties,
                            input_data.exp_conditions,comp_elements,comp_pst,input_data.element_data,
                            input_data.uncertainty_levels,input_data.structure)

        elif self.m.mode == "edge_single_calculation":
            print('---------Start Calculation for Selected Compositions----------')
            print('------------------Use Edge Dislocation Model------------------')
            input_data = read_inputjson_edge_single_calculation(self.inputfilename)
            input_data.check_integrity_curtin_edge()
            input_data.grab_properties_curtin_edge()
            self.input_data = input_data
            if input_data.uncertainty_levels == [0,0]:
                self.model = ss_edge_model_T(input_data.dislocation_properties,input_data.exp_conditions,
                            input_data.element_composition,input_data.element_data,input_data.structure)
            else: 
                self.model = ss_edge_model_T_w_uncertainty(ss_edge_model_T,input_data.dislocation_properties,
                            input_data.exp_conditions,input_data.element_composition,
                            input_data.element_data,input_data.uncertainty_levels,input_data.structure)

        elif self.m.mode == "screw_ternary":
            print('-------------Start Calculation for Ternary--------------')
            print('------------Use BCC-Screw Dislocation Model ------------')
            input_data = read_inputjson_BCC_screw_pseudo_ternary(self.inputfilename)
            self.input_data = input_data 
            mesh = build_mesh_ternary(input_data.increment,input_data.psA_range,input_data.psB_range,input_data.psC_range,
                          input_data.psA_elements,input_data.psB_elements,input_data.psC_elements,
                          input_data.psA_ratio,input_data.psB_ratio,input_data.psC_ratio)
            comp_tmp = mesh.make_mesh()
            comp_pst, comp_elements = mesh.comp_assign()
            self.model = ss_model_M_C_screw_pseudo_ternary(input_data,comp_elements,comp_pst)

        elif self.m.mode == "screw_single_calculation":
            print('---------Start Calculation for Selected Compositions----------')
            print('---------------Use BCC-Screw Dislocation Model----------------')
            input_data = read_inputjson_BCC_screw_single_calculation(self.inputfilename)
            self.input_data = input_data
            self.model = ss_model_M_C_screw(input_data)

        elif self.m.mode ==  "screw_suzuki_single_calculation":
            print('---------Start Calculation for Selected Compositions----------')
            print('------------Use Suzuki-BCC-Screw Dislocation Model------------')
            input_data = read_json_Suzuki_model_RWASM_T(self.inputfilename)
            self.input_data = input_data
            self.model = Suzuki_model_RWASM_T(input_data.element_data,input_data.experiment_conditions,
                                         input_data.adjustable_scalers)
        elif self.m.mode ==  "screw_suzuki_ternary_calculation":
            print('-------------Start Calculation for Ternary--------------')
            print('---------Use Suzuki-BCC-Screw Dislocation Model---------')
            input_data = read_json_Suzuki_model_RWASM_ternary(self.inputfilename)
            self.input_data = input_data
            mesh = build_mesh_ternary(input_data.increment,input_data.psA_range,input_data.psB_range,input_data.psC_range,
                          input_data.psA_elements,input_data.psB_elements,input_data.psC_elements,
                          input_data.psA_ratio,input_data.psB_ratio,input_data.psC_ratio)
            comp_tmp = mesh.make_mesh()
            comp_pst, comp_elements = mesh.comp_assign()
            self.model = Suzuki_model_RWASM_ternary(input_data.element_data,comp_elements,comp_pst,
                          input_data.experiment_conditions,input_data.adjustable_scalers)
        else: 
            print('invalid inputs. Please use valid format. check --help and --format')

        # start calculation and save data
        self.outputfile = self.input_data.savefilename
        print('-------------Running--------------')
        self.model.calculate()
        print('---------------Done----------------')
        self.model.writedata()
        print('-------------------------------Predicted Data-------------------------------')
        print(self.model.calc_data.reset_index(drop=True))
    def writeoutput(self):
        self.model.calc_data.to_csv(self.outputfile,index=False,sep='\t',float_format='%.2f')
        print('-------------------------------Data Saved to {}-------------------------------'.format(self.outputfile))

        if self.m.mode == "screw_suzuki_single_calculation":
            # write additional data to check convergence
            print('-------------------------------Additional Data for Suzuki Model-------------------------------')
            print('--------------------------Use this to check yield stress convergence--------------------------')
            for element_symbol in self.model.elements_kappa_i_convergence_record:
                print('Saving data for {}'.format(element_symbol))
                write_data_element = pd.DataFrame(data=self.model.elements_kappa_i_convergence_record[element_symbol].to_dict())
                write_data_element.to_csv(element_symbol+'_kappa_tau_y',index=False,sep='\t',float_format='%.5f')

