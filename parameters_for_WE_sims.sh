#!/bin/sh
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WE_string_exclusive_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.001,0.03,'exclusive_test_p_set_III_090215','III',1,4),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WE_string_exclusive_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.005,0.03,'exclusive_test_p_set_II_090215','II',1,4),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WE_string_exclusive_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.006,0.03,'exclusive_test_p_set_I_090215','I',1,4),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_090315('tri_gene_trial_1.net','WExplore_trigene_090315_trial_1.mat'),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_090315('tri_gene_trial_9.net','WExplore_trigene_090515_trial_9.mat'),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_090315('tri_gene_trial_10.net','WExplore_trigene_090515_trial_10.mat'),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_090315('tri_gene_activation_alt_7.net','WExplore_trigene_090515_alt_7.mat'),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_090315('tri_gene_activation_alt_6.net','WExplore_trigene_090515_alt_6.mat'),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_090315('tri_gene_activation_alt_7.net','WExplore_trigene_090815_alt_7_2.mat'),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_090315('tri_gene_activation_basal_3.net','WExplore_trigene_090815_basal_3.mat'),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WE_string_general_toggle_switch_ver_1(20,0.005,0.005,150,3000,10,5,0.006,0.03,'general_test_p_set_I_090215','I',1,4),catch,exit(1),end,exit(0)"

#/home/grad/BioNetGen-2.2.2-stable/bin/run_network_x86_64-linux -o ./tri_gene_activation_basal_4_long -p ssa -h $RANDOM --cdat 1 --fdat 0 -g ./tri_gene_activation_basal_4.net ./tri_gene_activation_basal_4.net 0.005 10000000
#/home/grad/BioNetGen-2.2.2-stable/bin/run_network_x86_64-linux -o ./tri_gene_activation_basal_4_long -p ssa -h $RANDOM --cdat 1 --fdat 0 -g ./tri_gene_activation_basal_5.net ./tri_gene_activation_basal_5.net 0.005 10000000

#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_091015('tri_gene_activation_basal_4.net','WExplore_trigene_091015_basal_4.mat'),catch,exit(1),end,exit(0)"
#/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_091015('tri_gene_activation_basal_5.net','WExplore_trigene_091015_basal_5_asymmetric.mat'),catch,exit(1),end,exit(0)"

/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_091015('tri_gene_activation_basal_6.net','WExplore_trigene_092115_basal_6.mat'),catch,exit(1),end,exit(0)"
/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_091015('tri_gene_activation_basal_7.net','WExplore_trigene_092115_basal_7_asymmetric.mat'),catch,exit(1),end,exit(0)"
/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_091015('tri_gene_activation_basal_8.net','WExplore_trigene_092115_basal_8.mat'),catch,exit(1),end,exit(0)"
/usr/local/MATLAB/R2012a/bin/matlab -nodisplay -r "try, WExplore_trigene_091015('tri_gene_activation_basal_9.net','WExplore_trigene_092115_basal_9_asymmetric.mat'),catch,exit(1),end,exit(0)"