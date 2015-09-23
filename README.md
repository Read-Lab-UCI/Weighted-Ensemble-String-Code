# WE-string-GRN

This repository contains the code used to generate the data presented in:
Tse et al., DNA-Binding Kinetics Determines the Mechanism of Noise-Induced Switching in Gene Networks, Biophysical
Journal (2015), http://dx.doi.org/10.1016/j.bpj.2015.08.035

The algorithm used in these files, the weighted ensemble-based string method was designed by:
Adelman, J. L., & Grabe, M. (2013). Simulating rare events using a weighted ensemble-based string method. 
The Journal of chemical physics, 138(4), 044105. http://dx.doi.org/10.1063/1.4773892

MATLAB version 2012b and BioNetGen ver 2.2.2 were used. The MATLAB parallel computing toolbox is required to run the parallelized portions of this code.

The parameters used to generate the data are as follows:
WE_string_exclusive_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.001,0.03,'output_filename1','III',1,4)
WE_string_exclusive_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.005,0.03,'output_filename2','II',1,4)
WE_string_exclusive_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.006,0.03,'output_filename3','I',1,4)
WE_string_general_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.001,0.03,'output_filename4','III',1,4)
WE_string_general_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.005,0.03,'output_filename5','II',1,4)
WE_string_general_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.006,0.03,'output_filename6','I',1,4)

One initial .bngl BioNetGen file: exclusive_switch_I.bngl is included.

The compiled BioNetGen files used in this paper are:
exclusive_switch_I.net	
exclusive_switch_II.net	
exclusive_switch_III.net	
general_switch_I.net	
general_switch_II.net	
general_switch_III.net

