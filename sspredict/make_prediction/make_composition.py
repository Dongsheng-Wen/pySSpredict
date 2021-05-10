import numpy as np
import pandas as pd 


class build_mesh_ternary:
##"build mesh grid for pseudo-ternary diagram"
# input is increment, and composition range for the pseudo ternary components
# generate evenly spaced concentration points for the pseudo ternary components (psA, psB, and psC)
# then for each pseudo ternary component, assign the elemental concentration according to the ratio

# e.g. psA = CoNi,
#      psA_elements: Co, Ni, 
#      psA_ratio: 1:1
#
    def __init__(self,increment,
                 psA_range,
                 psB_range,
                 psC_range,
                 psA_elements,
                 psB_elements,
                 psC_elements,
                 psA_ratio,
                 psB_ratio,
                 psC_ratio):
        if float(increment) <=0:
            print('error: increment should be positive float')
            quit()
        self.increment = float(increment)
        self.psA_range = psA_range
        self.psB_range = psB_range
        self.psC_range = psC_range
        
        self.psA_elements = psA_elements
        self.psB_elements = psB_elements
        self.psC_elements = psC_elements
        
        self.psA_ratio = np.array(psA_ratio)
        self.psB_ratio = np.array(psB_ratio)
        self.psC_ratio = np.array(psC_ratio)
        
        
    def make_mesh(self):
        self.comp_psA = []
        self.comp_psB = []
        self.comp_psC = []
        a1 = np.arange(0,100,self.increment);
        a1 = np.append(a1,100)
        for i in a1:
            a2 = np.arange(0,100-i,self.increment)
            a2 = np.append(a2,100-i)
            for j in a2:
                self.comp_psA.append(i)
                self.comp_psB.append(j)
                self.comp_psC.append(100-i-j)
        comp_tmp = pd.DataFrame(
            data={
                    "comp(psA)":self.comp_psA,
                    "comp(psB)":self.comp_psB,
                    "comp(psC)":self.comp_psC
                 })
        
        self.comp = comp_tmp[(comp_tmp['comp(psA)']>=self.psA_range[0]) & (comp_tmp['comp(psA)']<=self.psA_range[1])
                            &(comp_tmp['comp(psB)']>=self.psB_range[0]) & (comp_tmp['comp(psB)']<=self.psB_range[1])
                            &(comp_tmp['comp(psC)']>=self.psC_range[0]) & (comp_tmp['comp(psC)']<=self.psC_range[1])]
        
        return self.comp
    
    def comp_assign(self):
    # assign psA to psA_elements according to psA_ratio
    # assign psB to psB_elements according to psB_ratio
    # assign psC to psC_elements according to psC_ratio
    # e.g. 
    # if psA_elements = ['Co','Ni','Cr'], psA_ratio = [1,1,2]
    # when psA = 40 at.%, Co = 10 at.%, Ni = 10 at.%, Cr = 20 at.%
        self.comp_elements = pd.DataFrame(data={})
        self.comp_psA_elements = [self.comp['comp(psA)']*self.psA_ratio[i]/sum(self.psA_ratio) for i in range(len(self.psA_ratio))]
        self.comp_psB_elements = [self.comp['comp(psB)']*self.psB_ratio[i]/sum(self.psB_ratio) for i in range(len(self.psB_ratio))]
        self.comp_psC_elements = [self.comp['comp(psC)']*self.psC_ratio[i]/sum(self.psC_ratio) for i in range(len(self.psC_ratio))]
        
        for i in range(len(self.psA_elements)):
            self.comp[self.psA_elements[i]] = self.comp_psA_elements[i]
            self.comp_elements[self.psA_elements[i]] = self.comp_psA_elements[i]
        for i in range(len(self.psB_elements)):
            self.comp[self.psB_elements[i]] = self.comp_psB_elements[i]
            self.comp_elements[self.psB_elements[i]] = self.comp_psB_elements[i]
        for i in range(len(self.psC_elements)):
            self.comp[self.psC_elements[i]] = self.comp_psC_elements[i]
            self.comp_elements[self.psC_elements[i]] = self.comp_psC_elements[i]
        
        self.comp_new = self.comp.reset_index(drop=True)
        self.comp_elements_new = self.comp_elements.reset_index(drop=True)
        return self.comp_new, self.comp_elements_new