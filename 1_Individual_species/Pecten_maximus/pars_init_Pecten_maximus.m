function [par, metaPar, txtPar] = pars_init_Pecten_maximus(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 7671;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 1.05;       free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 12;         free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.030215;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.8865;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 49.06;    free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 2372*2;  free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 0.0002541; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 0.02141; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 3234; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 9.121e-09;  free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.del_M_P = 0.334;  free.del_M_P = 1;   units.del_M_P = '-';      label.del_M_P = 'shape coefficient for adult'; 
par.del_M_larv_P = 1.0667;  free.del_M_larv_P = 1;   units.del_M_larv_P = '-';  label.del_M_larv_P = 'shape coefficient for larvae'; 
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.f_tL = 0.4;   free.f_tL  = 1;   units.f_tL = '-';         label.f_tL = 'sc func res Bay of Brest, EVECOS'; 
par.f_tL2 = 0.6;        free.f_tL2 = 1;   units.f_tL2 = '-';        label.f_tL2 = 'sc func res bay of brest'; 
% par.f_tL1 = 1;        free.f_tL1 = 1;   units.f_tL1 = '-';        label.f_tL1 = 'sc func res Celtic sea'; 
% par.f_tL3 = 1;        free.f_tL3 = 1;   units.f_tL3 = '-';        label.f_tL3 = 'sc func res Norway'; 


%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 
par.d_V = 0.2; free.d_V = 0;   units.d_V = 'g dw/cm^3';        label.d_V = 'specific density of structure';
par.d_E = 0.2; free.d_E = 0;   units.d_E = 'g dw/cm^3';        label.d_E = 'specific density of reserve';
%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
