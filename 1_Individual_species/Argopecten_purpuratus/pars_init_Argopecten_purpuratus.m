function [par, metaPar, txtPar] = pars_init_Argopecten_purpuratus(metaData)

metaPar.model = 'abj'; 

%% core primary parameters 
par.z = 0.9754;       free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 12;         free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.03876;      free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.9719;        free.kap   = 0;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
% par.kap = 0.8;        free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 51.68;      free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 2374*2;       free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 2.859e-06;  free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 0.0004535;  free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 75.57;     free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 5.544e-08;  free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.T_A = 8129;         free.T_A   = 0;         units.T_A = 'K';            label.T_A = 'Arrhenius temperature'; 
par.T_ref = 293.15;     free.T_ref = 0;         units.T_ref = 'K';          label.T_ref = 'Reference temperature'; 
par.del_M_A = 0.3743;	free.del_M_A = 0;       units.del_M_A = '-';        label.del_M_A = 'shape coefficient Argopecten'; 
par.del_M_larv_A = 0.7;	free.del_M_larv_A = 0;  units.del_M_larv_A = '-';	label.del_M_larv_A = 'shape coefficient for larvae Argopecten'; 
par.f = 1;              free.f = 0;             units.f = '-';              label.f = 'scaled functional response for 0-var data'; 
par.f_Paracas = 0.95; 	free.f_Paracas = 1;     units.f_Paracas = '-';      label.f_Paracas = 'sc func res Paracas Bay'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);
par.d_V = 0.2; free.d_V = 0;   units.d_V = 'g dw/cm^3';        label.d_V = 'specific density of structure';
par.d_E = 0.2; free.d_E = 0;   units.d_E = 'g dw/cm^3';        label.d_E = 'specific density of reserve';
%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
