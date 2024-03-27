function [par, metaPar, txtPar] = pars_init_group(metaData)
global pets % add global pets to be able to get the right values in metaData
metaPar.model = {'abj', 'abj', 'abj', 'abj', 'abj'}; 
% order of species:
% 'Argopecten_purpuratus', 'Mimachlamys_varia', 'Nodipecten_subnodosus', 'Pecten_maximus', 'Placopecten_magellanicus'
metaPar.cov_rules = ''; % see function parGrp2Pets
metaPar.weights.p_M = 0; % 1 % 2 % 3 % 4 % 5 % 6
metaPar.weights.v = 0; % 1 % 2 % 3 % 4 % 5 % 6


%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% calculation for body size scaling
%%% before adding the body size scaling
%par.del_M_A = 0.28022;  free.del_M_A = 1;   units.del_M_A = '-';      label.del_M_A = 'shape coefficient Argopecten purpuratus'; 
%par.del_M_M = 0.25031;  free.del_M_M = 1;   units.del_M_M = '-';      label.del_M_M = 'shape coefficient Mimachlamys varia'; 
%par.del_M_N = 0.30616;  free.del_M_N = 1;   units.del_M_N = '-';      label.del_M_N = 'shape coefficient Nodpiecten subnodosus'; 
%par.del_M_P = 0.28182;  free.del_M_P = 1;   units.del_M_P = '-';      label.del_M_P = 'shape coefficient Pecten maximus'; 
%par.del_M_Pl = 0.3014;  free.del_M_Pl = 1;   units.del_M_Pl = '-';     label.del_M_Pl = 'shape coefficient Placopecten magellanicus'; 

% hypothesis of all shape coefficient equal (mean at 0.28 from above)
par.del_M_A = 0.28;  free.del_M_A = 0;   units.del_M_A = '-';      label.del_M_A = 'shape coefficient Argopecten purpuratus'; 
par.del_M_M = 0.28;  free.del_M_M = 0;   units.del_M_M = '-';      label.del_M_M = 'shape coefficient Mimachlamys varia'; 
par.del_M_N = 0.28;  free.del_M_N = 0;   units.del_M_N = '-';      label.del_M_N = 'shape coefficient Nodpiecten subnodosus'; 
par.del_M_P = 0.28;  free.del_M_P = 0;   units.del_M_P = '-';      label.del_M_P = 'shape coefficient Pecten maximus'; 
par.del_M_Pl = 0.28;  free.del_M_Pl = 0;   units.del_M_Pl = '-';     label.del_M_Pl = 'shape coefficient Placopecten magellanicus'; 

% creation of vector of ultimate length for each species, from ultimate length in data saved in metaData to have access to the value
par.Li_species = [metaData.(pets{1}).Li metaData.(pets{2}).Li metaData.(pets{3}).Li metaData.(pets{4}).Li metaData.(pets{5}).Li];  free.Li_species   = [0 0 0 0 0];   units.Li_species = 'cm';   label.Li_species = 'ultimate shell height per species'; 
free.Li_species   = [0 0 0 0 0];   
units.Li_species = 'cm';   label.Li_species = 'ultimate shell height per species'; 

% zoom factor given for the initial set of parameters. Free, and estimated during the procedure
par.z_ref = 0.5489; free.z_ref     = 1;   units.z_ref = '-';            label.z_ref = 'zoom factor';
% calculation of the "species zoom factor" to represent the ratio between ultimate size of species and the ultimate size of the reference species
par.z_species = [(par.Li_species(1)*par.del_M_A)/(par.Li_species(5) *   par.del_M_Pl) ... 
                (par.Li_species(2) * par.del_M_M)/(par.Li_species(5) * par.del_M_Pl) ...
                (par.Li_species(3) * par.del_M_N)/(par.Li_species(5) * par.del_M_Pl) ...
                (par.Li_species(4) * par.del_M_P)/(par.Li_species(5) * par.del_M_Pl) ... 
                (par.Li_species(5) * par.del_M_Pl)/(par.Li_species(5) * par.del_M_Pl)]; % here 5 is P. magellanicus, the reference species, thus z_species = 1 for Pl
free.z_species = [0 0 0 0 0]; units.z_species = '-'; label.z_species = 'species zoom factor, ratio between ultimate sizes';


%% core primary parameters 
par.z = par.z_ref * par.z_species; free.z= [0 0 0 0 0]; units.z = '-'; label.z = 'initial zoom factor per species';
par.F_m = 6.5;          free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;        free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;        free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v =  [0.0203 0.0203 0.0203 0.0203 0.0203]; %0.0203; % % depends if we want to have different values or not
free.v     =  [1 1 1 1 1]; %1; %
units.v = 'cm/d';         label.v = 'energy conductance';
par.kap =  0.9232; free.kap   = 1;  units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;       free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M =     [35.0041 35.0041 35.0041 35.0041 35.0041];      %35.0041; % % depends if we want to have different values or not
free.p_M   =  [1 1 1 1 1];   %1;%
units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;            free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;        free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 4185;         free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 3.1052e-4;    free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hj = 0.3248;        free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp_ref = 1123;   free.E_Hp_ref  = 1;   units.E_Hp_ref = 'J';         label.E_Hp_ref = 'reference maturity at puberty for body size scaling'; 
par.E_Hp = [ par.E_Hp_ref*((par.z_species(1)).^3) par.E_Hp_ref*((par.z_species(2)).^3) par.E_Hp_ref*((par.z_species(3)).^3) par.E_Hp_ref*((par.z_species(4)).^3) par.E_Hp_ref*((par.z_species(5)).^3)]; free.E_Hp  = [0 0 0 0 0 ];   units.E_Hp = 'J';         label.E_Hp = 'initial maturity at puberty with body size scaling'; 
        free.E_Hp  = [0 0 0 0 0];   units.E_Hp = 'J';         label.E_Hp = 'initial maturity at puberty with body size scaling'; 
par.h_a = 5.0524e-09;     free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;       free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

par.T_A = 8078;         free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 

%% other parameters 
par.T_AH = 4168.539;        free.T_AH  = 0;   units.T_AH = 'K';     label.T_AH = 'Arrhenius temperature above the upper boundary'; 
par.T_AL = 176544.7305;     free.T_AL  = 0;   units.T_AL = 'K';     label.T_AL = 'Arrhenius temperature below the lower boundary'; 
par.T_H = 293.15;           free.T_H   = 0;   units.T_H = 'K';      label.T_H = 'upper boundary tolerance range'; 
par.T_L = 275.15;           free.T_L   = 0;   units.T_L = 'K';      label.T_L = 'lower boundary tolerance range'; 

% par.del_M_A = 0.28022;  free.del_M_A = 1;   units.del_M_A = '-';      label.del_M_A = 'shape coefficient Argopecten purpuratus'; 
% par.del_M_M = 0.25031;  free.del_M_M = 1;   units.del_M_M = '-';      label.del_M_M = 'shape coefficient Mimachlamys varia'; 
% par.del_M_N = 0.30616;  free.del_M_N = 1;   units.del_M_N = '-';      label.del_M_N = 'shape coefficient Nodpiecten subnodosus'; 
% par.del_M_P = 0.28182;  free.del_M_P = 1;   units.del_M_P = '-';      label.del_M_P = 'shape coefficient Pecten maximus'; 
% par.del_M_Pl = 0.3014;  free.del_M_Pl = 1;   units.del_M_Pl = '-';     label.del_M_Pl = 'shape coefficient Placopecten magellanicus'; 

par.del_M_larv_A = 0.5212;  free.del_M_larv_A = 0;   units.del_M_larv_A = '-';  label.del_M_larv_A = 'shape coefficient for larvae Argopecten purpuratus'; 
par.del_M_larv_M = 0.8296;  free.del_M_larv_M = 0;   units.del_M_larv_M = '-';  label.del_M_larv_M = 'shape coefficient for larvae Mimachlamys varia'; 
par.del_M_larv_N = 1.043;   free.del_M_larv_N = 0;   units.del_M_larv_N = '-';  label.del_M_larv_N = 'shape coefficient for larvae Nodipecten subnodosus'; 
par.del_M_larv_P = 1.07;    free.del_M_larv_P = 0;   units.del_M_larv_P = '-';  label.del_M_larv_P = 'shape coefficient for larvae Pecten maximus'; 
par.del_M_larv_Pl = 0.7942; free.del_M_larv_Pl = 0;  units.del_M_larv_Pl = '-'; label.del_M_larv_Pl = 'shape coefficient for larvae Placopecten magellanicus'; 
par.f = 1;              free.f     = 0;     units.f = '-';            label.f =         'scaled functional response for 0-var data'; 
par.f_Paracas = 1;      free.f_Paracas = 1; units.f_Paracas = '-';    label.f_Paracas = 'sc func res Paracas Bay (A. purpuratus, Aguirre-Velarde et al., 2019)'; 
par.f_tL19SA =  1.1971;       free.f_tL19SA = 1;  units.f_tL19SA = '-';     label.f_tL19SA =  'sc func res Bay of Brest (M. varia, Régnier-Brisson et al. in prep)'; 
par.f_AM =      1.0063;         free.f_AM  = 1;     units.f_AM = '-';         label.f_AM =      'sc func res laguna Ojo de liebre (N.subnodosus, Arellano-Martinez 2011)'; 
par.f_CL =      1.1770;      free.f_CL  = 1;     units.f_CL = '-';         label.f_CL =      'sc func res laguna Ojo de liebre and Bahia de Los Angeles (N. subnodosus, Carreno-Leon et al., 2023)'; 
par.f_RA =      0.5406;        free.f_RA  = 1;     units.f_RA = '-';         label.f_RA =      'sc func res Bahia de Loreto (N. subnodosus, Ramirez-Arce 2009)'; 
par.f_RC =      0.9836;         free.f_RC  = 1;     units.f_RC = '-';         label.f_RC =      'sc func res Bahia Magdalena (N. subnodosus, Racotta et al., 2003)'; 
par.f_MA =      0.8345;        free.f_MA  = 1;     units.f_MA = '-';         label.f_MA =      'sc func res Bahia Magdalena, high density (N. subnodosus, Maldonado-Amparo et al., 2004)'; 
par.f_EV =      1.0358;          free.f_EV  = 1;     units.f_EV = '-';         label.f_EV =      'sc func res Bay of Brest (P. maximus, EVECOS)'; 
par.f_B =       1.1427;         free.f_B = 1;     units.f_B = '-';        label.f_B =     'sc func res Bay of Brest (P. maximus, Lavaud 2014)'; 
par.f_CB =      0.9642;         free.f_CB = 1;     units.f_CB = '-';        label.f_CB =     'sc func res Baie des Chaleurs (P. magellanicus, Claereboudt et al., 1994)'; 
par.f_MD10 =    1.1756;       free.f_MD10 = 1;   units.f_MD10 = '-';      label.f_MD10 =   'sc func res Newfoundland, 10 m deep (P. magellanicus, Macdonald & Thompson 1988)'; 
par.f_MD31 = 	0.94;       free.f_MD31 = 1;   units.f_MD31 = '-';      label.f_MD31 =   'sc func res Newfoundland, 31 m deep (P. magellanicus, Macdonald & Thompson 1988)'; 
par.f_RO =      0.79;        free.f_RO = 1;    units.f_RO = '-';       label.f_RO =    'sc func res Gulf of Saint-Lawrence (P. magellanicus, Roddick et al., 1994)'; 

%% set chemical parameters from Kooy2010 
[phylum, class] = metaData2taxo(metaData); 
[par, units, label, free] = addchem(par, units, label, free, phylum, class); 
% modify the specific density of structure to consider less water content
par.d_V = 0.2; free.d_V = 0;   units.d_V = 'g dw/cm^3';        label.d_V = 'specific density of structure';
par.d_E = 0.2; free.d_E = 0;   units.d_E = 'g dw/cm^3';        label.d_E = 'specific density of reserve';

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
