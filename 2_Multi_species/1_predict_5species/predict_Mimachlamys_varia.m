function [prdData, info] = predict_Mimachlamys_varia(par, data, auxData)
  
  % unpack par, data, auxData
    %%% calculation for body size scaling relationships
    par.E_Hp = par.E_Hp_ref*((par.z_species)^3);
    cPar = parscomp_st(par);
    cPar.p_Am = (par.z_ref * par.z_species) * par.p_M/ par.kap;   % J/d.cm^2, {p_Am} spec assimilation flux; the expression for p_Am is multiplied also by L_m^ref = 1 cm, for units to match. 

  vars_pull(par); % unpacks variables from the structure where the primary parameters are
  vars_pull(cPar);  % unpacks variables from the structure where the compound parameters are
  vars_pull(data); % unpacks variables from the structure where the data are
  vars_pull(auxData); % unpacks variables from the structure where the auxiliary data are
  
  % compute temperature correction factors
  % for each data with temperature given
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_aj = tempcorr(temp.aj, T_ref, T_A);
  TC_tp = tempcorr(temp.tp, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  %TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  TC_RL = tempcorr(temp.R_L, T_ref, T_A);
%   TC_19SA = tempcorr(temp.tL19SA, T_ref, T_A); 
%   TC_LW19SA = tempcorr(temp.LW19SA, T_ref, T_A);
%   TC_tLlarvae = tempcorr(temp.tLlarvae, T_ref, T_A);
   
   filterChecks = ...% f contrained to not be larger than 1.2 or negative
               f_tL19SA > 1.2 || f_tL19SA < 0 ;% || ...; % 
%                del_M_N < del_M_P || del_M_P < del_M_A || del_M_A < del_M_Pl || del_M_Pl < del_M_M;
if filterChecks  
    info = 0;
    prdData = {};
    return;
end

%%
    % -----------------------------------------------------------------------
    % zero-variate data
    % -----------------------------------------------------------------------  
  % as zero-variate data for Lp Li are issued from L. Regnier-Brisson's in
  % situ monitoring, we will use f_tL19SA for f. But in fact use f = 1 for
  % all during the multi-species
  
   % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v];                % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0);       % d.cm^2, initial scaled reserve
  E_0 = J_E_Am * U_E0 * mu_E;                       % J, initial reserve for kap and v constant
 
%    f = f_Tinduff;% Hatchery for early life stages
%   f = f_tL19SA;

  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
  
  % birth
  L_b = L_m * l_b;                          % cm, structural length at birth
  Lw_b = L_b/ del_M_larv_A;                 % cm, shell height at birth
  aT_b = t_b/ k_M/ TC_ab;               	% d, age at birth

  % metamorphosis
  L_j = L_m * l_j;                          % cm, structural length at metam
  Lw_j = L_j/ del_M_M;                      % cm, shell height at metam
  tT_j = t_j / k_M/ TC_aj;                  % d, time since fertilisation at metam 
  
%    f = f_tL19SA;% 2019 Sainte Anne
%   [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

  % puberty 
  L_p = L_m * l_p;                          % cm, structural length at puberty
  Lw_p = L_p/ del_M_M;                      % cm, shell height at puberty
  tT_p = (t_p - t_b)/ k_M/ TC_tp;           % d, time since birth at puberty
  Ww_p = L_p.^3  *(1 + f * w );             % g, wet weight
  
  % ultimate
  L_i = L_m * l_i;                          % cm, ultimate structural length 
  Lw_i = L_i/ del_M_M  ;                    % cm, ultimate total length 
  Wd_i = d_V * L_i.^3  *(1 + f * w );       % g, dry weight without gonad

 
  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];      % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);          % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/TC_am;                    % d, mean life span at T
  
   % reproduction (from P. maximus)
%   f = f_repro;% reproduction conditions (mix between Tinduff hatchery and Saint-Anne monitoring)
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp];     % compose parameter vector at T
  LR = L_R.R_L;                                                     % call the right size to calculate reproductive rate
  RT_L = TC_RL * reprod_rate_j(LR, f, pars_R);                      % #/d, ultimate reproduction rate at T

  E_R = RT_L * 365/2 * E_0 / kap_R;                                 % J, recalculate reproduction buffer from reproductive rate
  W_at_LR = d_V * (LR * del_M_M).^3  *(1 + f * w );                 % g, dry weight without gonad at f
  W_ER_at_LR = w_E / mu_E * E_R ;                                   % g, weight of reproduction buffer at reproduction size given
  GSI = W_ER_at_LR / (W_at_LR + W_ER_at_LR);                        % _, gonadosomatic index
 

  % pack to output
  prdData.ab = aT_b;    % d, age at birth
  prdData.aj = tT_j;    % d, age at metamorphosis since fertilisation
  prdData.tp = tT_p;    % d, time since birth at puberty at f and T     
  prdData.am = aT_m;    % d, life span
  prdData.Lb = Lw_b;    % cm, length at birth
  prdData.Lj = Lw_j;    % cm, length at metamorphosis
  prdData.Lp = Lw_p;    % cm, shell height at puberty at f              
  prdData.Li = Lw_i;    % cm, ultimate length
  prdData.Wwp = Ww_p;   % g, wet weight at puberty
  prdData.Wdi = Wd_i;   % g, ultimate dry weight
%   prdData.Ri = RT_i;    % # d^-1, ultimate reproduction rate
  prdData.R_L = RT_L;    % # d^-1, reproduction rate at L_R
  prdData.GSI = GSI;    % %, gonado-somatic index
  prdData.E_0 = E_0;    % J, energy on an egg

 %%
  % -----------------------------------------------------------------------
  % uni-variate data
  % -----------------------------------------------------------------------
  
  %%% Time-Length of larvae and juveniles, after birth and before puberty
  %------------------------------------------------------------------------
  % Data from Tinduff hatchery
%   [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_Tinduff); %f_Tinduff, the functional response for all data sets coming from Tinduff hatchery
%   
%   %%% Larvae 
%   rT_j = TC_tLlarvae * rho_j * k_M;         % exponential growth rate, corrected with the temperature
%   time_bj = tLlarvae(:,1);                  % days, time between birth and metamorphosis
%   Lw_b = l_b * L_m/ del_M_larv_A;           % cm, physical length at birth
%   ELw_bj = Lw_b * exp(time_bj * rT_j/3);    % estimated length between birth and metamorphosis
%   
%   %%% Juveniles
%   rT_B = TC_tLlarvae * rho_B * k_M;         % von Bertalanffy growth rate, corrected with the temperature
%   time_ji = tLjuve(:,1);                    % days, time after metamorphosis
%   aT_j = t_j/ k_M/ TC_tLlarvae;             % time at metamorphosis, since fertilisation
%   Lw_j = l_j * L_m/ del_M_M;                % cm, physical length at metamorphosis
%   Lw_i = l_i * L_m/ del_M_M;                % cm, physical ultimate length
%   ELw_ji = Lw_i - (Lw_i - Lw_j) * exp( - rT_B * (time_ji - aT_j));  % estimated length after metamorphosis
  
  
  %------------------------------------------------------------------------
  %%% Time-Length-Ash free dry weight, of juveniles and adults
  %------------------------------------------------------------------------
  % Saint-Anne/Roscanvel data set, from Laure Regnier-Brisson thesis
  % Temperature values from environmental data set, defined here with
  % Fourier transform function and parameters.
 
  [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_tL19SA);
  L_b = l_b * L_m;   L_j = l_j * L_m;   % cm, structural length at birth and at metamorphosis
  s_M = L_j/ L_b;                       % -, acceleration factor at f  
  
  tT = temp.tLh_19SA;                     % K, temperature from data
  L0 = Lw0.tLh_19SA * del_M_M;            % cm, initial structural length defined in data
  f = f_tL19SA;                         % -, scaled function response used to compute initial condition and to get the dynamic
  E0 = f * E_m * L0^3;                  % J, initial amount of energy in reserves
  VE0  = [L0^3; E0];                    % [cm^3; J], initial values of state variables (volume and reserve)
  
  [t_sort_19SA, it , it_sort_19SA] = unique(tLh_19SA(:,1),'sorted'); %take the time of the data
  
  [~, VE] =  ode45(@dget_VE, t_sort_19SA, VE0,[], f, tT, par, cPar, s_M);       % resolution of the dynamic with differential equations
  Lw   = (VE(:,1).^(1/3)) / del_M_M;                                            % cm, physical shell height
  ELw_19SA = Lw(it_sort_19SA);                                                  % cm, reconstruction, physical shell height for data time
  
  % to consider time of weight data set for custom plots
  [t_sort_19SA, it , it_sort_19SA] = unique(tWd_19SA(:,1),'sorted');
  [~, VE] =  ode45(@dget_VE, t_sort_19SA, VE0,[], f, tT, par, cPar, s_M);
  W_s = d_V * (VE(:,1).^(1/3)).^3 + VE(:,2) .* w_E / mu_E / d_E;                 % g, somatic dry weight
  EtWd_19SA = W_s(it_sort_19SA);                                                   % g dw, reconstruction, somatic dry weight for data time

  %length-weight estimation /!\ E_R not considered
  EWd_19SA = d_V * ((LWd_19SA(:,1)*del_M_M ).^3  *(1 + f_tL19SA * w )) ;             %in g, dry weight of tissues and not include reproduction buffer

    
    
  %------------------------------------------------------------------------
  % pack to output
%   prdData.tLlarvae = ELw_bj;
%   prdData.tLjuve = ELw_ji;
  prdData.tLh_19SA = ELw_19SA;  
  prdData.tWd_19SA = EtWd_19SA;
  prdData.LWd_19SA = EWd_19SA;

end
  
%%% subfunctions
%--------------------------------------------------------------------------
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
        pS = p.p_M * V + p.p_T * V^(2/3);                                       % J d^-1, somatic maintenance rate
        pG = p.kap * pC - pS;                                                   % J d^-1, growth rate

        % Differential equations of state variables  
        dE = pA - pC;               % J d^-1, reserve dynamic
        dV = pG / p.E_G;            % J d^-1, volume dynamic

        % Pack output 
        dVE = [dV; dE]; 
end
     