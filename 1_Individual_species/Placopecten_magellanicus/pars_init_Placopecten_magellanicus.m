function [par, metaPar, txtPar] = pars_init_Placopecten_magellanicus(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 9000;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 0.5046;       free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.016847;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.913;      free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 31.204;     free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5215;       free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 3.011e-04; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 5.924e-01; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metamorphosis'; 
par.E_Hp = 1.026e+03; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 4.321e-09;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.T_AH = 4168.539;  free.T_AH  = 0;   units.T_AH = 'K';         label.T_AH = 'Arrhenius temperature above the upper boundary'; 
par.T_AL = 176544.7305;  free.T_AL  = 0;   units.T_AL = 'K';         label.T_AL = 'Arrhenius temperature below the lower boundary'; 
par.T_H = 293.15;     free.T_H   = 0;   units.T_H = 'K';          label.T_H = 'upper boundary tolerance range'; 
par.T_L = 275.15;     free.T_L   = 0;   units.T_L = 'K';          label.T_L = 'lower boundary tolerance range'; 
par.del_M_Pl = 0.30535;  free.del_M_Pl = 1;   units.del_M_Pl = '-';     label.del_M_Pl = 'shape coefficient'; 
par.del_Mb = 0.80412;  free.del_Mb = 1;   units.del_Mb = '-';       label.del_Mb = 'shape coefficient before metamorphosis'; 
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.f_LWd = 1.2;      free.f_LWd = 1;   units.f_LWd = '-';        label.f_LWd = 'scaled functional response for LWd data'; 
par.f_tL1 = 1.1386;   free.f_tL1 = 1;   units.f_tL1 = '-';        label.f_tL1 = 'scaled functional response for tL data at 10 m deep in Newfoundland'; 
par.f_tL2 = 1.1989;   free.f_tL2 = 1;   units.f_tL2 = '-';        label.f_tL2 = 'scaled functional response for tL data at 31 m deep in Newfoundland'; 
par.f_tL3 = 0.75701;  free.f_tL3 = 1;   units.f_tL3 = '-';        label.f_tL3 = 'scaled functional response for tL data in the Gulf of Saint-Lawrence'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 
par.d_V = 0.2;        free.d_V   = 0;   units.d_V = 'g dw/cm^3';  label.d_V = 'specific density of structure'; 
par.d_E = 0.2;        free.d_E   = 0;   units.d_E = 'g dw/cm^3';  label.d_E = 'specific density of reserve'; 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
