function derivative = flux(t, EVHR, simu)

%------------------------------------------------------
% Objectives: calculate the fluxes at each time and the derivative of the
% state variable. It is called by ode function to realise the integration
% Only for adult stage

% Inputs:   t for the time
%           EVHR for the initial values of state variable
%           simu all information for each individual
%
% Outputs: derivatives of E, V, EH and ER
%
% calls:        food.m
%               temperature.m
%
% called by:    indiv.m

% 2024/02/06 - Eline Le Moan - based on Laure Pecquerie's codes
% Mascoet Project
%------------------------------------------------------

vars_pull(simu.par)
vars_pull(simu.parC)

% state variables
E = EVHR(1);            % J,    reserve
V = EVHR(2);            % cm^3, structure
E_H = EVHR(3);          % J,    maturity
E_R = EVHR(4);          % J,    reproduction buffer


%% Environmental conditions
% calculation of forcing variables from temperature and food value
f = food(t, simu) / (food(t, simu) + simu.X_K);
cT = exp((T_A / T_ref)-(T_A/temperature(t, simu)));

% correction for temperature
v   =   v * cT;
p_M =   p_M * cT;
k_J =   k_J * cT;
p_Am =  p_Am * cT;

% % metabolic acceleration
p_Am    = p_Am * simu.parStat.s_M;
v       = v * simu.parStat.s_M;

% body size
p_Am = p_Am * simu.z_species;
E_Hp = E_Hp * simu.z_species^3;

%% Flux calculations
% assimilation flux
if E_H < E_Hb
    p_A = 0;
else
    p_A = p_Am * f * V^(2/3);
end
% somatic maintenance flux
p_S = p_M * V;
% mobilisation flux
p_C = (E/V) * (E_G * v * V^(2/3) + p_S) / (kap * (E/V) + E_G);
% growth flux
p_G = kap * p_C - p_S;
% maturity maintenance flux
p_J = k_J * E_H;
% reproduction flux
p_R = (1 - kap) * p_C - p_J;

%% Derivative calculations
% variation of reserve
dE = p_A - p_C;
% variation of structure
dV = p_G / E_G;
% variation of maturity and reproduction buffer
if E_H < E_Hp
    dEH = p_R;
    dER = 0;
elseif E_H >= E_Hp 
    dEH = 0;
    dER = kap_R * p_R;
end

derivative = [dE; dV; dEH; dER]; % derivative value for each state variable
end