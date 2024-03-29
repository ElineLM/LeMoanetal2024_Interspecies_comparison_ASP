function [par, metaPar, txtPar] = pars_init_Nodipecten_subnodosus(metaData)

metaPar.model = 'abj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 6577;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 1.258;       free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 12;         free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.01479;     free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.9857;    free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 19.32;    free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 5308;     free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 1.33e-05; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 0.002572; free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 1681; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 1.76e-08;  free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.del_M_N = 0.3710;  free.del_M_N = 1;   units.del_M_N = '-';      label.del_M_N = 'shape coefficient'; 
par.del_M_larv_N = 1.0346;  free.del_M_larv_N = 0;   units.del_M_larv_N = '-';  label.del_M_larv_N = 'shape coefficient for larvae'; 
par.f = 1;              free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.f_AM = 1.2;         free.f_AM  = 1;   units.f_AM = '-';         label.f_AM = 'sc func res laguna Ojo de liebre (Arellano-Martinez 2011)'; 
par.f_CL = 0.4869;      free.f_CL  = 1;   units.f_CL = '-';         label.f_CL = 'sc func res laguna Ojo de liebre and Bahia de Los Angeles (Carreno-Leon et al., 2023)'; 
par.f_RA = 0.7372;       free.f_RA  = 1;   units.f_RA = '-';         label.f_RA = 'sc func res Bahia de Loreto (Ramirez-Arce 2009)'; 

par.f_RC = 1.2;     free.f_RC  = 1;   units.f_RC = '-';         label.f_RC = 'sc func res  (Racotta et al., 2003)'; 
par.f_MA = 0.9811;     free.f_MA  = 1;   units.f_MA = '-';         label.f_MA = 'sc func res  (Maldonado-Amparo et al., 2004)'; 
% par.f_NG = 0.6;     free.f_NG  = 1;   units.f_NG = '-';         label.f_NG = 'sc func res  ()'; 
% par.f_SG = 0.6;     free.f_SG  = 1;   units.f_SG = '-';         label.f_SG = 'sc func res  ()'; 
%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 
par.d_V = 0.2;        free.d_V   = 0;   units.d_V = 'g dw/cm^3';  label.d_V = 'specific density of structure'; 
par.d_E = 0.2;        free.d_E   = 0;   units.d_E = 'g dw/cm^3';  label.d_E = 'specific density of reserve'; 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
