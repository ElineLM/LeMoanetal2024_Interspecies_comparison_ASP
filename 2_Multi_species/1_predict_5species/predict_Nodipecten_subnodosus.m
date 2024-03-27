function [prdData, info] = predict_Nodipecten_subnodosus(par, data, auxData)
  % unpack par, data, auxData
  %%% calculation for body size scaling relationships
    par.E_Hp = par.E_Hp_ref*((par.z_species)^3);
    cPar = parscomp_st(par);
    cPar.p_Am = (par.z_ref * par.z_species) * par.p_M/ par.kap;   % J/d.cm^2, {p_Am} spec assimilation flux; the expression for p_Am is multiplied also by L_m^ref = 1 cm, for units to match. 
    vars_pull(par);     vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
 % del_M max fixed at 1.3 as the shape coeff for a sphere
 filterChecks = ...% f contrained to not be larger than 1.2 or negative
               f_RA > 1.2 || f_RA < 0 || ...;
               f_AM > 1.2 || f_AM < 0 || ...;
               f_CL > 1.2 || f_CL < 0 || ...;
               f_RC > 1.2 || f_RC < 0 || ...;
               f_MA > 1.2 || f_MA < 0 ;% || ...;
%                del_M_N < del_M_P || del_M_P < del_M_A || del_M_A < del_M_Pl || del_M_Pl < del_M_M;
if filterChecks  
    info = 0;
    prdData = {};
    return;
end  
  
  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_tj = tempcorr(temp.tj, T_ref, T_A);
  TC_tp = tempcorr(temp.tp, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
%   TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  TC_RL = tempcorr(temp.R_L, T_ref, T_A);
%   TC_tLh_RC = tempcorr(temp.tLh_RC, T_ref, T_A);
%   TC_NG22  = tempcorr(temp.NG_tLl, T_ref, T_A);
%   TC_RA    = tempcorr(C2K(25), T_ref, T_A);
%   TC_GWw   = tempcorr(temp.tGWw, T_ref, T_A);

%%
% -----------------------------------------------------------------------  
% zero-variate data
% -----------------------------------------------------------------------
    % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

   % to calculate E_0, considering same f than reproduction
  pars_UE0 = [V_Hb; g; k_J; k_M; v];            % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0);   % d.cm^2, initial scaled reserve
  E_0 = J_E_Am * U_E0 * mu_E;                   % J, initial reserve for kap and v constant

  % birth
  L_b  = L_m * l_b;                             % cm, structural length at birth at f
  Lw_b = L_b / del_M_larv_N;                    % cm, shell height at birth at f
  Wd_b = L_b^3 * d_V * (1 + f * w);             % g, dry weight at birth
  aT_b = t_b / k_M / TC_ab;                     % d, age at birth at f and T

  % metam
  L_j  = L_m * l_j;                             % cm, structural length at metam
  Lw_j = L_j/ del_M_N;                          % cm, shell height at metam at f
  Wd_j = L_j^3 * d_V * (1 + f * w);             % g, dry weight at metam
  tT_j = t_j / k_M / TC_tj;                     % d, time since birth at metam
  s_M  = l_j / l_b;                             % - , acceleration factor before metamorphosis
  
  % puberty 
  L_p = L_m * l_p;                              % cm, structural length at puberty at f
  Lw_p = L_p/ del_M_N;                          % cm, shell height at puberty at f
  Ww_p = L_p^3 * (1 + f * w);                   % g, wet weight at puberty 
  tT_p = (t_p - t_b)/ k_M/ TC_tp;               % d, time since birth at puberty at f and T

  % ultimate
  L_i = L_m * l_i;                              % cm, ultimate structural length at f
  Lw_i = L_i/ del_M_N;                          % cm, ultimate shell height at f
  Ww_i = L_i^3 * (1 + f * w);                   % g, ultimate dry weight 
 
  % reproduction
%   pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
%   RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % #/d, ultimate reproduction rate at T
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
  LR = L_R.R_L;
  RT_L = TC_RL * reprod_rate_j(LR, f, pars_R);                 % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  % pack to output
  prdData.ab = aT_b;
  prdData.tj = tT_j;
  prdData.tp = tT_p;
  prdData.am = aT_m;
  prdData.Lb = Lw_b;
  prdData.Lj = Lw_j;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Wdb = Wd_b;
  prdData.Wdj = Wd_j;
  prdData.Wwp = Ww_p;
  prdData.Wwi = Ww_i;
%   prdData.Ri = RT_i;
  prdData.R_L = RT_L;
  prdData.E_0 = E_0;
  
%%
  % -----------------------------------------------------------------------
  % uni-variate data
  % -----------------------------------------------------------------------
  
  % Racotta 2003, time(d)~length(cm) 
%   t = tLh_RC(:,1); % d, time since birth
%   f = f_RC; % -, mean scaled func response for Racotta 2003
%   [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
%   L_b = l_b * L_m; 
%   L_j = l_j * L_m;
%   L_i = l_i * L_m; % cm, struct. length at birth
%   kT_M = k_M * TC_tLh_RC;
%   rT_j = rho_j * kT_M; 
%   rT_B = rho_B * kT_M; % 1/d, rates
%   tT_j = (t_j - t_b)/ kT_M; % d, time end of acceleration
%   L_bj = L_b * exp(t(t < tT_j,1) * rT_j/ 3);
%   L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t(t >= tT_j) - tT_j)); % cm, expected length at time
%   ELw_RC = [L_bj; L_jm]/ del_M_N;                                       
  %%% modified on March 2024 to use ode function and specify initial values
    f = f_RC;                      % -, initial scaled function response is used to compute initial condition of the individual
    [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f);
    L_b = l_b * L_m;   L_j = l_j * L_m;                         % cm, structural length at birth and at metamorphosis
    s_M = L_j/ L_b;                                             % -, acceleration factor at f
    time = tLh_RC(:,1);                                         % d, time from data
    tT = temp.tLh_RC;                                           % K, temperature from data
    L0 = Lw0.tLh_RC * del_M_N;                                  % cm, initial structural length
    E0 = f * E_m * L0^3;                                        % J, initial amount of energy in reserves
    VE0 = [L0^3; E0];                                           % [cm^3, J], initial values of the state variables (structure and reserve)
    [~, VE] =  ode45(@dget_VE, time, VE0,[],f, tT, par, cPar, s_M); % dynamic of the state variables, structure (V) and reserve (E)
    ELw_RC   = (VE(:,1).^(1/3)) / del_M_N;                      % cm, physical shell height


  
%   % food conditions Ramirez-Arce 2009
%   x = food.tWw_RA(:,2)  / K_ram; % scaled amount of food in environment
%   tf = x ./ (x + 1);  % scaled functional response for Ramirez-Arce 2009
%   tf = [food.tWw_RA(:,1), tf]; % appended times to scaled food
%   
%   % Ramirez-Arce 2009, time(d)~ biomass weight(g) ~shell height (cm) 
    time = tLh_RA(:,1);                                       % d, time from data
    tT = temp.tLh_RA;                                         % K, temperature from data
    tf = f_RA;                                                % -, scaled functional response                
    L0 = Lw0.tLh_RA * del_M_N;                                % cm, initial structural length
    f_init = f_RA;                                            % -, initial scaled functional response
    E0 = f_init * E_m * L0^3;                                 % J, initial amount of energy in reserves
    VE0 = [L0^3; E0];                                         % [cm^3, J], initial values of the state variables (structure and reserve)
    sM = L_j/L_b;                                             % -, acceleration factor at f
    [~, VE] = ode45(@dget_VE, time, VE0, [], tf, tT, par, cPar, sM); % dynamic of the state variables, structure (V) and reserve (E)
    L = VE(:,1).^(1/3);                                       % cm, structural shell height
    ELw_RA = L / del_M_N;                                     % cm, physical shell height
    % if include wet weight at time
    % modified to have the time for weight data set for custom graphs
    time = tWw_RA(:,1);
    [~, VE] = ode45(@dget_VE, time, VE0, [], tf, tT, par, cPar, sM);
    E = VE(:,2);
    EWw_RA = VE(:,1) + E * w_E / mu_E;                        % g, expected biomass wet weight, made the hypothesis of low ER weight

%     % wet weight at shell height -- Ramirez-Arce 2009
%     ELWw_RA = (LWw_RA(:,1)*del_M_N ).^3  *(1 + f_RA * w ) ; %in g, wet weight

%   % Nava-Gomez 2022, time(d)~length(cm) for larvae and juveniles because
%   % from 0 to 45 and metamorphosis at 15 days
%   f = f_NG;
%   t = NG_tLl(:,1);
%   [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
%   L_b = l_b * L_m; 
%   L_j = l_j * L_m;
%   L_i = l_i * L_m; % cm, struct. length 
%   kT_M = k_M * TC_NG22;
%   rT_j = rho_j * kT_M; 
%   rT_B = rho_B * kT_M; % 1/d, rates
%   tT_j = (t_j - t_b)/ kT_M; % d, time end of acceleration
%   L_bj = L_b * exp(t(t < tT_j,1) * rT_j/ 3);
%   L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t(t >= tT_j) - tT_j)); % cm, expected length at time
%   ELw_NG = [L_bj; L_jm]/ del_M_N;    
%   
%   
%   % Serrano-Guzman 1997{2} -- only larvae because before the time of
%   % metamorphosis
%   f = f_SG;
%   [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
%   L_b  = L_m * l_b;  Lw_b_SG = L_b / del_M_larv_N;                 % cm, structural / physical length at birth at f_SG
%   rT_j =  rho_j * k_M;
%   ELw_SG = Lw_b_SG * exp((SG_tLl(:,1) - aT_b) * rT_j/ 3);  % cm, expected length of settled larva
  
  % Arellano-Martinez 2011
    time = tLh_AM(:,1);
    tT = temp.tLh_AM;
    tf = f_AM;
    L0 = Lw0.tLh_AM * del_M_N;
    f_init = f;
    E0 = f_init * E_m * L0^3;
    VE0 = [L0^3; E0];
    sM = L_j/L_b;
    [~, VE] = ode45(@dget_VE, time, VE0, [], tf, tT, par, cPar, sM);
    L = VE(:,1).^(1/3);
    ELw_AM = L / del_M_N;
  

   % Carreño-Leon et al., 2023
    EWd_CL = d_V * ((LWd_CL(:,1)*del_M_N ).^3  *(1 + f_CL * w )) ; %in g, dry weight, without ER
  
     % Maldonado-Amparo et al., 2004
    time = tLl_MA(:,1);
    tT = temp.tLl_MA;
    tf = f_MA;
    L0 = Lw0.tLl_MA * del_M_N;
    f_init = f;
    E0 = f_init * E_m * L0^3;
    VE0 = [L0^3; E0];
    sM = L_j/L_b;
    
    [t_sort_AM, it , it_sort_AM] = unique(tLl_MA(:,1),'sorted'); %take the time of the data
    [~, VE] = ode45(@dget_VE, t_sort_AM, VE0, [], tf, tT, par, cPar, sM);
    Lw = (VE(:,1).^(1/3)) / del_M_N;
    ELw_MA = Lw(it_sort_AM);
    
  % pack to output
  prdData.tLh_RC    = ELw_RC;
  prdData.tWw_RA    = EWw_RA;
%   prdData.LWw_RA    = ELWw_RA;
  prdData.tLh_RA    = ELw_RA; 
%   prdData.NG_tLl    = ELw_NG;
%   prdData.SG_tLl    = ELw_SG;
  prdData.tLh_AM    = ELw_AM;
   prdData.LWd_CL    = EWd_CL;
   prdData.tLl_MA    = ELw_MA;
  
  
end
  
%modified on 19/09/2023 - Eline Le Moan
function dVE = dget_VE(t, VE, f, tT, p, c, s_M) 
    % INPUTS: time, initial values (volume and reserve), food, temperature, primary parameters, compound parameters, acceleration factor
    % OUTPUTS: derivative of volume and reserve
    % modified from function in Pecten maximus
    
        V = VE(1);     % cm^3, structural volume
        E = VE(2);     % J, reserve
        
        if length(tT) > 1 % if temperature is a table
            if tT(2,1) < 10 && tT(2,2) < 10 % search if value of time and temperature are below 10, it means that this corresponds to fourier parameters
                tT_int = fnfourier(t, tT)+273.15;      
            else 
                tT_int = interp1(tT(:,1),tT(:,2),t,'linear');  % interpolate the right temperature
            end 
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
