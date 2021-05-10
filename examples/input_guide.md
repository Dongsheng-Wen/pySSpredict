line1: title line  
line2: increment of composition (0 to 100)  
line3: strain rate  
line4: temperature  
line5: - group elements into pseudo-ternary components (psA,psB,psC),  
       &ensp;&ensp;&ensp;&ensp;&ensp;- separated by ',', within one component separated by '-'. example 'e: Mn-Co, Fe-Ni,Al'  
line6: - set ratios for grouped elements in pseudo-ternary components (rA,rB,rC)  
 &ensp;&ensp;&ensp;&ensp;&ensp;- separated by ',', numbers are separated by '-' within one component,  
 &ensp;&ensp;&ensp;&ensp;&ensp;- if this component contains only one element: use '1  
 &ensp;&ensp;&ensp;&ensp;&ensp;- int, float will work. example 'ratio:1-1,1-3,1', 'ratio:0.5-0.4-0.1,10-90,1'  
line7: structure of the alloy, it can be either fcc or bcc  
line8: 'data:' indicate the start of material data  
line9: element names separated by ','  
line10: atomic volumes of elements, separated by ',', unit: Angstrom  
line11:Young's Modulus of elements, separated by ',', unit: GPa  
line12:Shear Modulus of elements, separated by ',', unit: GPa  

----
an example is provided below, more input examples can be found in examples/  
MnFeCoNiAl  
increment: 10  
strain_r: 0.001  
temperature: 300  
e: Mn-Co, Fe-Ni,Al  
ratio: 1-1,1-2,1  
structure: fcc  
data:  
Mn,Fe,Co,Ni,Al  
12.60,12.09,11.12,10.94,16.472  
338.84,26.375,262.89,199.1,65.498  
130.52,8.9429,101.66,76.0,23.922  
