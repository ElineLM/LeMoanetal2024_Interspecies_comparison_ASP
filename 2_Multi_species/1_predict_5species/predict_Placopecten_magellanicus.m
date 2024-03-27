function [prdData, info] = predict_Placopecten_magellanicus(par, data, auxData)

% Unpack par, data, auxData
%%% calculation for body size scaling relationship
    par.E_Hp = par.E_Hp_ref*((par.z_species)^3);
    cPar = parscomp_st(par);
    cPar.p_Am = (par.z_ref * par.z_species) * par.p_M/ par.kap;  % J/d.cm^2, {p_Am} spec assimilation flux; the expression for p_Am is multiplied also by L_m^ref = 1 cm, for units to match. 
    vars_pull(par); vars_pull(cPar); vars_pull(data); vars_pull(auxData);

filterChecks = ...%
    f_MD10 < 0 || f_MD10 > 1.2 || ...;
    f_MD31 < 0 || f_MD31 > 1.2 || ...;
    f_RO < 0  || f_RO > 1.2  || ...;
    f_CB < 0   || f_CB > 1.2     ;% ||...
%     del_M_N < del_M_P || del_M_P < del_M_A || del_M_A < del_M_Pl || del_M_Pl < del_M_M;

if filterChecks
  prdData = []; 
  info = 0; 
  return
end

%% Compute temperature correction factors
Tpars    = [T_A T_L T_H T_AL T_AH];
TC_ab    = tempcorr(temp.ab, T_ref, Tpars);
TC_aj    = tempcorr(temp.aj, T_ref, Tpars);
TC_ap    = tempcorr(temp.ap, T_ref, Tpars);
TC_am    = tempcorr(temp.am, T_ref, Tpars);
% TC_Ri    = tempcorr(temp.Ri, T_ref, Tpars);
TC_RL    = tempcorr(temp.R_L, T_ref, Tpars);
TC_MD10   = tempcorr(temp.tLh_MD10, T_ref, Tpars);
TC_MD31   = tempcorr(temp.tLh_MD31, T_ref, Tpars);
% TC_RO   = tempcorr(temp.tLh_RO, T_ref, Tpars);
% TC_LN    = tempcorr(temp.LN, T_ref, Tpars);
% TC_WJO   = tempcorr(temp.WJO, T_ref, Tpars);
% TC_WJO1  = tempcorr(temp.WJO1, T_ref, Tpars);
% TC_WJO3  = tempcorr(temp.WJO3, T_ref, Tpars);
% TC_WJO6  = tempcorr(temp.WJO6, T_ref, Tpars);
% TC_WJO8  = tempcorr(temp.WJO8, T_ref, Tpars);
% TC_WJO10 = tempcorr(temp.WJO10, T_ref, Tpars);
% TC_WJO19 = tempcorr(temp.WJO19, T_ref, Tpars);
% TC_TJO   = tempcorr(C2K(TJO(:,1)), T_ref, Tpars);

%% 
        % -----------------------------------------------------------------------
        % Zero-variate data
        % -----------------------------------------------------------------------
 % Life cycle
pars_tj = [g k l_T v_Hb v_Hj v_Hp];
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, ~, ~, info] = get_tj(pars_tj, f); % -, scaled times & lengths at f
if info == 0
  prdData = []; return;
end

    % to calculate E_0, considering same f than reproduction
    pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
    U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
    E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap and v constant

	% Birth
    L_b  = L_m * l_b;                 % cm, structural length at birth at f
    Lw_b = L_b / del_M_larv_Pl;             % cm, physical length at birth at f
    Wd_b = L_b^3 * d_V * (1 + f * w); % g,  dry weight at birth at f (remove d_V for wet weight)
    aT_b = t_b / k_M / TC_ab;         % d,  age at birth at f and T

    % Metamorphosis
    L_j  = L_m * l_j;                 % cm, structural length at metamorphosis at f
    Lw_j = L_j / del_M_larv_Pl;             % cm, physical length at metamorphosis at f
    Wd_j = L_j^3 * d_V * (1 + f * w); % g,  dry weight at metamorphosis at f (remove d_V for wet weight)
    aT_j = t_j / k_M / TC_aj;         % d,  age at metamorphosis at f and T

     % Puberty
    L_p  = L_m * l_p;                 % cm, structural length at puberty at f
    Lw_p = L_p / del_M_Pl;               % cm, physical length at puberty at f
    Wd_p = L_p^3 * d_V * (1 + f * w); % g,  dry weight at puberty (remove d_V for wet weight)
    aT_p = t_p / k_M / TC_ap;         % d,  age at puberty at f and T

     % Ultimate
    L_i  = L_m * l_i;                 % cm, ultimate structural length at f
    Lw_i = L_i / del_M_Pl;               % cm, ultimate physical length at f
    Wd_i = L_i^3 * d_V * (1 + f * w); % g,  ultimate dry weight (remove d_V for wet weight)

     % Reproduction
    % pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
    % RT_i   = TC_Ri * reprod_rate_j(L_i, f, pars_R);               % #/d, ultimate reproduction rate at T
    pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
    LR = L_R.R_L;
    RT_L = TC_RL * reprod_rate_j(LR, f, pars_R);                 % #/d, ultimate reproduction rate at T

     % Life span
    pars_tm = [g; l_T; h_a / k_M^2; s_G]; % compose parameter vector at T_ref
    t_m     = get_tm_s(pars_tm, f, l_b);  % -, scaled mean life span at T_ref
    aT_m    = t_m / k_M / TC_am;          % d, mean life span at T

     % Pack to output
    prdData.ab  = aT_b;
    prdData.aj  = aT_j;
    prdData.ap  = aT_p;
    prdData.am  = aT_m;
    prdData.Lb  = Lw_b;
    prdData.Lj  = Lw_j;
    prdData.Lp  = Lw_p;
    prdData.Li  = Lw_i;
    prdData.Wdb = Wd_b;
    prdData.Wdj = Wd_j;
    prdData.Wdp = Wd_p;
    prdData.Wdi = Wd_i;
    % prdData.Ri  = RT_i;
    prdData.R_L = RT_L;
    prdData.E_0 = E_0;

%% 
        % -----------------------------------------------------------------------
        % Uni-variate data
        % -----------------------------------------------------------------------
     % Time vs. shell length at 10 m deep in Newfoundland
%     f = f_MD10; TC = TC_MD10; tL = tLh_MD10;
%     [t_j, ~, t_b, l_j, ~, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
%     rT_j = TC * rho_j * k_M; rT_B = TC * rho_B * k_M;
%     L_b  = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;
%     tT_j = (t_j - t_b) / k_M / TC;                                   % d,   time since birth at metamorphosis
%     t_bj = tL(tL(:,1) < tT_j,1);                                     % d,   select times before metamorphosis
%     L_bj = (L_b * exp(t_bj * rT_j / 3)) / del_M_larv_Pl;                    % cm,  physical length at exponential growth as V1-morph
%     t_ji = tL(tL(:,1) >= tT_j,1);                                    % d,   select times after metamorphosis
%     L_ji = (L_i - (L_i - L_j) * exp(-rT_B * (t_ji - tT_j))) / del_M_Pl; % cm,  physical length at isomorphic growth
%     EL1  = [L_bj; L_ji];                                             % cm,  expected physical length at time
    %%% modified on March 2024 to use ode function and specify initial values
    [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_MD10);
      L_b = l_b * L_m;   L_j = l_j * L_m; % cm, structural length at birth and at metamorphosis
     s_M = L_j/ L_b;                 % -, acceleration factor at f
    tT = temp.tLh_MD10;                  % K, temperature from data
    L0 = Lw0.tLh_MD10 * del_M_Pl;           % cm, initial structural length
    f = f_MD10;                      % -, initial scaled function response is used to compute initial condition of the individual
    E0 = f * E_m * L0^3;            % J, initial amount of energy in reserves
    VE0 = [L0^3; E0];               % [cm^3, J], initial values of the state variables (structure and reserve)
    [t_sort_MD10, it , it_sort_MD10] = unique(tLh_MD10(:,1),'sorted'); %take the time of the data
    [~, VE] =  ode45(@dget_VE, t_sort_MD10, VE0,[],f, tT, par, cPar, s_M);
    Lw   = (VE(:,1).^(1/3)) / del_M_Pl;           % cm, physical shell height
    EL1 = Lw(it_sort_MD10);     % cm, physical shell height with all time
    
    % 31 m deep in Newfoundland
%     f = f_MD31; TC = TC_MD31; tL = tLh_MD31;
%     [t_j, ~, t_b, l_j, ~, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
%     rT_j = TC * rho_j * k_M; rT_B = TC * rho_B * k_M;
%     L_b  = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;
%     tT_j = (t_j - t_b) / k_M / TC;                                   % d,   time since birth at metamorphosis
%     t_bj = tL(tL(:,1) < tT_j,1);                                     % d,   select times before metamorphosis
%     L_bj = (L_b * exp(t_bj * rT_j / 3)) / del_M_larv_Pl;                    % cm,  physical length at exponential growth as V1-morph
%     t_ji = tL(tL(:,1) >= tT_j,1);                                    % d,   select times after metamorphosis
%     L_ji = (L_i - (L_i - L_j) * exp(-rT_B * (t_ji - tT_j))) / del_M_Pl; % cm,  physical length at isomorphic growth
%     EL2  = [L_bj; L_ji];                                             % cm,  expected physical length at time
    %%% modified on March 2024 to use ode function and specify initial values
    [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_MD31);
      L_b = l_b * L_m;   L_j = l_j * L_m; % cm, structural length at birth and at metamorphosis
     s_M = L_j/ L_b;                 % -, acceleration factor at f
    tT = temp.tLh_MD31;                  % K, temperature from data
    L0 = Lw0.tLh_MD31 * del_M_Pl;           % cm, initial structural length
    f = f_MD31;                      % -, initial scaled function response is used to compute initial condition of the individual
    E0 = f * E_m * L0^3;            % J, initial amount of energy in reserves
    VE0 = [L0^3; E0];               % [cm^3, J], initial values of the state variables (structure and reserve)
    [t_sort_MD31, it , it_sort_MD31] = unique(tLh_MD31(:,1),'sorted'); %take the time of the data
    [~, VE] =  ode45(@dget_VE, t_sort_MD31, VE0,[],f, tT, par, cPar, s_M);
    Lw   = (VE(:,1).^(1/3)) / del_M_Pl;           % cm, physical shell height
    EL2 = Lw(it_sort_MD31);     % cm, physical shell height with all time
 
    
    % 31 m deep in New Jersey
    % f = f_tL3; TC = TC_tL3; tL = tL3;
    % [t_j, ~, t_b, l_j, ~, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
    % rT_j = TC * rho_j * k_M; rT_B = TC * rho_B * k_M;
    % L_b  = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;
    % tT_j = (t_j - t_b) / k_M / TC;                                   % d,   time since birth at metamorphosis
    % t_bj = tL(tL(:,1) < tT_j,1);                                     % d,   select times before metamorphosis
    % L_bj = (L_b * exp(t_bj * rT_j / 3)) / del_M_larv_Pl;                    % cm,  physical length at exponential growth as V1-morph
    % t_ji = tL(tL(:,1) >= tT_j,1);                                    % d,   select times after metamorphosis
    % L_ji = (L_i - (L_i - L_j) * exp(-rT_B * (t_ji - tT_j))) / del_M_Pl; % cm,  physical length at isomorphic growth
    % EL3  = [L_bj; L_ji];                                             % cm,  expected physical length at time

    % Time vs. shell length at 10 m deep in the Gulf of Saint-Lawrence
    % f = f_tLSL; TC = TC_tLSL; tL = tLSL;
    % [t_j, ~, t_b, l_j, ~, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
    % rT_j = TC * rho_j * k_M; rT_B = TC * rho_B * k_M;
    % L_b  = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;
    % tT_j = (t_j - t_b) / k_M / TC;                                   % d,   time since birth at metamorphosis
    % t_bj = tL(tL(:,1) < tT_j,1);                                     % d,   select times before metamorphosis
    % L_bj = (L_b * exp(t_bj * rT_j / 3)) / del_M_larv_Pl;                    % cm,  physical length at exponential growth as V1-morph
    % t_ji = tL(tL(:,1) >= tT_j,1);                                    % d,   select times after metamorphosis
    % L_ji = (L_i - (L_i - L_j) * exp(-rT_B * (t_ji - tT_j))) / del_M_Pl; % cm,  physical length at isomorphic growth
    % EL3  = [L_bj; L_ji];                                             % cm,  expected physical length at time
    %%% modified on March 2024 to use ode function and specify initial values
    [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_RO);
      L_b = l_b * L_m;   L_j = l_j * L_m; % cm, structural length at birth and at metamorphosis
     s_M = L_j/ L_b;                 % -, acceleration factor at f
    tT = temp.tLh_RO;                  % K, temperature from data
    L0 = Lw0.tLh_RO * del_M_Pl;           % cm, initial structural length
    f = f_RO;                      % -, initial scaled function response is used to compute initial condition of the individual
    E0 = f * E_m * L0^3;            % J, initial amount of energy in reserves
    VE0 = [L0^3; E0];               % [cm^3, J], initial values of the state variables (structure and reserve)
    [t_sort_SL, it , it_sort_SL] = unique(tLh_RO(:,1),'sorted'); %take the time of the data
    [~, VE] =  ode45(@dget_VE, t_sort_SL, VE0,[],f, tT, par, cPar, s_M);
    Lw   = (VE(:,1).^(1/3)) / del_M_Pl;           % cm, physical shell height
    ELwSL = Lw(it_sort_SL);     % cm, physical shell height with all time


    % Shell length vs. tissue dry weight
    f   = f_CB;
    EWd = (LWd_CB(:,1) * del_M_Pl).^3 * (1 + w * f) * d_V; % g, expected dry weight
    % do not consider reproduction buffer in it

    %  % Shell length vs. fecundity
    % f  = f_LN; TC = TC_LN;
    % EN = 365 * TC * reprod_rate_j(LN(:,1) * del_M_Pl, f, pars_R); % #, number of offspring per year
    % 
    % % Compute parameters for oxygen consumption rate predictions
    % p_ref  = p_Am * L_m^2; % J/d, max assimilation power at max size and T_ref
    % pars_p = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp];
    % 
    %  % Dry weight vs. oxygen consumption rate
    % f = f_WJO; TC = TC_WJO; W = WJO(:,1);
    % [~, ~, ~, l_j, l_p, l_b, ~, ~, ~] = get_tj(pars_tj, f);
    % L        = (W / d_V / (1 + f * w)).^(1/3);
    % pACSJGRD = TC * p_ref * scaled_power_j(L, f, pars_p, l_b, l_j, l_p);
    % pADG     = pACSJGRD(:,[1 7 5])';   % J/d,    assimilation, dissipation, growth power
    % J_O      = eta_O * pADG;           % mol/d,  fluxes of organics J_X, J_V, J_E, J_P in rows
    % J_M      = -n_M \ n_O * J_O;       % mol/d,  fluxes of minerals J_C, J_H, J_O, J_N in rows, A, D, G in cols
    % EWJO     = -32e3 / 24 * J_M(3,:)'; % mgO2/h, O2 consumption
    % 
    %  % Dry weight vs. oxygen consumption rate at 1 °C
    % f = f_WJO1; TC = TC_WJO1; W = WJO1(:,1);
    % [~, ~, ~, l_j, l_p, l_b, ~, ~, ~] = get_tj(pars_tj, f);
    % L        = (W / d_V / (1 + f * w)).^(1/3);
    % pACSJGRD = TC * p_ref * scaled_power_j(L, f, pars_p, l_b, l_j, l_p);
    % pADG     = pACSJGRD(:,[1 7 5])';     % J/d,    assimilation, dissipation, growth power
    % J_O      = eta_O * pADG;             % mol/d,  fluxes of organics J_X, J_V, J_E, J_P in rows
    % J_M      = -n_M \ n_O * J_O;         % mol/d,  fluxes of minerals J_C, J_H, J_O, J_N in rows, A, D, G in cols
    % EJO1     = -24.4e3 / 24 * J_M(3,:)'; % mLO2/h, O2 consumption
    %  % 3 °C
    % f = f_WJO3; TC = TC_WJO3; W = WJO3(:,1);
    % [~, ~, ~, l_j, l_p, l_b, ~, ~, ~] = get_tj(pars_tj, f);
    % L        = (W / d_V / (1 + f * w)).^(1/3);
    % pACSJGRD = TC * p_ref * scaled_power_j(L, f, pars_p, l_b, l_j, l_p);
    % pADG     = pACSJGRD(:,[1 7 5])';     % J/d,    assimilation, dissipation, growth power
    % J_O      = eta_O * pADG;             % mol/d,  fluxes of organics J_X, J_V, J_E, J_P in rows
    % J_M      = -n_M \ n_O * J_O;         % mol/d,  fluxes of minerals J_C, J_H, J_O, J_N in rows, A, D, G in cols
    % EJO3     = -24.4e3 / 24 * J_M(3,:)'; % mLO2/h, O2-consumption
    %  % 6 °C
    % f = f_WJO6; TC = TC_WJO6; W = WJO6(:,1);
    % [~, ~, ~, l_j, l_p, l_b, ~, ~, ~] = get_tj(pars_tj, f);
    % L        = (W / d_V / (1 + f * w)).^(1/3);
    % pACSJGRD = TC * p_ref * scaled_power_j(L, f, pars_p, l_b, l_j, l_p);
    % pADG     = pACSJGRD(:,[1 7 5])';     % J/d,    assimilation, dissipation, growth power
    % J_O      = eta_O * pADG;             % mol/d,  fluxes of organics J_X, J_V, J_E, J_P in rows
    % J_M      = -n_M \ n_O * J_O;         % mol/d,  fluxes of minerals J_C, J_H, J_O, J_N in rows, A, D, G in cols
    % EJO6     = -24.4e3 / 24 * J_M(3,:)'; % mLO2/h, O2-consumption
    %  % 8 °C
    % f = f_WJO8; TC = TC_WJO8; W = WJO8(:,1);
    % [~, ~, ~, l_j, l_p, l_b, ~, ~, ~] = get_tj(pars_tj, f);
    % L        = (W / d_V / (1 + f * w)).^(1/3);
    % pACSJGRD = TC * p_ref * scaled_power_j(L, f, pars_p, l_b, l_j, l_p);
    % pADG     = pACSJGRD(:,[1 7 5])';     % J/d,    assimilation, dissipation, growth power
    % J_O      = eta_O * pADG;             % mol/d,  fluxes of organics J_X, J_V, J_E, J_P in rows
    % J_M      = -n_M \ n_O * J_O;         % mol/d,  fluxes of minerals J_C, J_H, J_O, J_N in rows, A, D, G in cols
    % EJO8     = -24.4e3 / 24 * J_M(3,:)'; % mLO2/h, O2-consumption
    %  % 10 °C
    % f = f_WJO10; TC = TC_WJO10; W = WJO10(:,1);
    % [~, ~, ~, l_j, l_p, l_b, ~, ~, ~] = get_tj(pars_tj, f);
    % L        = (W / d_V / (1 + f * w)).^(1/3);
    % pACSJGRD = TC * p_ref * scaled_power_j(L, f, pars_p, l_b, l_j, l_p);
    % pADG     = pACSJGRD(:,[1 7 5])';     % J/d,    assimilation, dissipation, growth power
    % J_O      = eta_O * pADG;             % mol/d,  fluxes of organics J_X, J_V, J_E, J_P in rows
    % J_M      = -n_M \ n_O * J_O;         % mol/d,  fluxes of minerals J_C, J_H, J_O, J_N in rows, A, D, G in cols
    % EJO10    = -24.4e3 / 24 * J_M(3,:)'; % mLO2/h, O2-consumption
     % 19 °C
    % f = f_WJO19; TC = TC_WJO19; W = WJO19(:,1);
    % [~, ~, ~, l_j, l_p, l_b, ~, ~, ~] = get_tj(pars_tj, f);
    % L        = (W / d_V / (1 + f * w)).^(1/3);
    % pACSJGRD = TC * p_ref * scaled_power_j(L, f, pars_p, l_b, l_j, l_p);
    % pADG     = pACSJGRD(:,[1 7 5])';     % J/d,    assimilation, dissipation, growth power
    % J_O      = eta_O * pADG;             % mol/d,  fluxes of organics J_X, J_V, J_E, J_P in rows
    % J_M      = -n_M \ n_O * J_O;         % mol/d,  fluxes of minerals J_C, J_H, J_O, J_N in rows, A, D, G in cols
    % EJO19    = -24.4e3 / 24 * J_M(3,:)'; % mLO2/h, O2-consumption
    % % 
    %  % Temperature vs. oxygen consumption rate
    % f = f_TJO; W = Wd0;
    % [~, ~, ~, l_j, l_p, l_b, ~, ~, ~] = get_tj(pars_tj, f);
    % L        = (W / d_V / (1 + f * w)).^(1/3);
    % pACSJGRD = p_ref * scaled_power_j(L, f, pars_p, l_b, l_j, l_p);
    % pADG     = pACSJGRD(:,[1 7 5])';                     % J/d,    assimilation, dissipation, growth power
    % J_O      = eta_O * pADG;                             % mol/d,  fluxes of organics J_X, J_V, J_E, J_P in rows
    % J_M      = -n_M \ n_O * J_O;                         % mol/d,  fluxes of minerals J_C, J_H, J_O, J_N in rows, A, D, G in cols
    % ETJO     = -24.4e3 / 24 * J_M(3,:)' .* TC_TJO / Wd0; % mLO2/h, O2-consumption

     % Pack to output
    prdData.tLh_MD10   = EL1;
    prdData.tLh_MD31   = EL2;
    prdData.tLh_RO   = ELwSL;
    prdData.LWd_CB   = EWd;
    % prdData.LN    = EN;
    % prdData.WJO   = EWJO;
    % prdData.WJO1  = EJO1;
    % prdData.WJO3  = EJO3;
    % prdData.WJO6  = EJO6;
    % prdData.WJO8  = EJO8;
    % prdData.WJO10 = EJO10;
    % prdData.WJO19 = EJO19;
    % prdData.TJO   = ETJO;

end


%modified on 19/09/2023 - Eline Le Moan
function dVE = dget_VE(t, VE, f, tT, p, c, s_M) 
    % INPUTS: time, initial values (volume and reserve), food, temperature, primary parameters, compound parameters, acceleration factor
    % OUTPUTS: derivative of volume and reserve
    % modified from function in Pecten maximus
    
        V = VE(1);     % cm^3, structural volume
        E = VE(2);     % J, reserve
        
        if length(tT) > 1
            tT_int = interp1(tT(:,1),tT(:,2),t,'linear');  % interpolate the right temperature
        else
            tT_int = tT; % get the temperature value
        end
        TC = tempcorr(tT_int, p.T_ref, p.T_A); % -, temperature correction factor

        % Shape correction function applies to surface-area specific assimilation and energy conductance
        % all parameters with time in dimension need to be corrected by the temperature
        p.v    = s_M * p.v * TC;            % cm d^-1, conductance corrected
        c.p_Am = s_M * c.p_Am * TC;         % J cm^-2 d^-1, maximum specific assimilation rate corrected
        p.p_M  = p.p_M * TC;                % J cm^-3 d^-1, volume-specific somatic maintenance rate

        % Rates
        pA  = f * c.p_Am * V^(2/3) ;                                            % J d^-1, assimilation rate
        pC  = E * (p.E_G * p.v/ V^(1/3)  + p.p_M)/ (p.kap * E/V + p.E_G);       % J d^-1, mobilisation rate
        pS = p.p_M * V + p.p_T * V^(2/3);                                      % J d^-1, somatic maintenance rate
        pG = p.kap * pC - pS;                                                % J d^-1, growth rate

        % Differential equations of state variables  
        dE = pA - pC;               % J d^-1, reserve dynamic
        dV = pG / p.E_G;           % J d^-1, volume dynamic

        % Pack output 
        dVE = [dV; dE]; 
end