function [prdData, info] = predict_Argopecten_purpuratus(par, data, auxData)
  
    % unpack par, data, auxData
    
    %%% calculation for body size scaling relationships
%    par.E_Hp = par.E_Hp*((par.z_species)^3);
    cPar = parscomp_st(par);
%   cPar.p_Am = (par.z_ref * par.z_species) * par.p_M/ par.kap;   % J/d.cm^2, {p_Am} spec assimilation flux; the expression for p_Am is multiplied also by L_m^ref = 1 cm, for units to match. 
    
    vars_pull(par); vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);


    filterChecks = ...% f contrained to not be negative and above 1.2
               f_Paracas < 0 || f_Paracas > 1.2 ;% || ...
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
%     TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
    TC_RL = tempcorr(temp.R_L, T_ref, T_A);
%     TC_tL = tempcorr(temp.tL, T_ref, T_A);
%     TC_tLlarv = tempcorr(temp.tLlarv, T_ref, T_A) ;

%%
    % -----------------------------------------------------------------------
    % zero-variate data
    % -----------------------------------------------------------------------
    % life cycle
    pars_tj = [g k l_T v_Hb v_Hj v_Hp];
    [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

     % to calculate E_0, considering same f than reproduction
    pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
    U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
    E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap and v constant

    % birth
    L_b = L_m * l_b;                            % cm, structural length at birth
    Lw_b = L_b/ del_M_larv_A;                 % cm, shell height at birth
    Wd_b = L_b^3 * d_V * (1 + f * w);           % g, dry weight at birth
    aT_b = t_b/ k_M/ TC_ab;                     % d, age at birth

    % metam
    L_j = L_m * l_j;                            % cm, structural length at metam
    Lw_j = L_j/ del_M_A;                        % cm, shell height at metam
    Wd_j = L_j^3 * d_V * (1 + f * w);           % g, dry weight at metam
    tT_j = (t_j - t_b)/ k_M/ TC_tj;             % d, time since birth at metam

    % puberty 
    L_p = L_m * l_p;                            % cm, structural length at puberty
    Lw_p = L_p/ del_M_A;                        % cm, shell height at puberty
    Wd_p = L_p^3 * d_V * (1 + f * w);           % g, dry weight at puberty 
    tT_p = (t_p - t_b)/ k_M/ TC_tp;             % d, time since birth at puberty

    % ultimate
    L_i = L_m * l_i;                            % cm, ultimate structural length
    Lw_i = L_i/ del_M_A;                        % cm, ultimate shell height
    Wd_i = L_i^3 * d_V * (1 + f * w);           % g, ultimate dry weight (without reproduction buffer)

    % reproduction
%     pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
%     RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % #/d, ultimate reproduction rate at T
    pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T
    LR = L_R.R_L;
    RT_L = TC_RL * reprod_rate_j(LR, f, pars_R);                 % #/d, ultimate reproduction rate at T

    % life span
    pars_tm = [g; l_T; h_a/ k_M^2; s_G];        % compose parameter vector at T_ref
    t_m = get_tm_s(pars_tm, f, l_b);            % -, scaled mean life span at T_ref
    aT_m = t_m/ k_M/ TC_am;                     % d, mean life span at T

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
    prdData.Wdp = Wd_p;
    prdData.Wdi = Wd_i;
%     prdData.Ri = RT_i;
    prdData.R_L = RT_L;
    prdData.E_0 = E_0;
    
%%   
    % -----------------------------------------------------------------------
    % uni-variate data
    % -----------------------------------------------------------------------
    
%     % time-length univariate data tL after metamorphosis (> 80 days)
%     % ------------------------------
%     %     f = f_tL;             % if A. purpuratus alone
%     %     f = f_Paracas;        % if multispecies estimation
%     [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f_Paracas); % directly written f_Paracas
%     t = tL(:,1);                                                            % d, time with birth at zero, from data
%     kT_M = k_M * TC_tL; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;           % 1/d, rates
%     tT_j = (t_j - t_b)/ kT_M;                                               % d, time of the end of acceleration
%     L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t - tT_j));                   % cm, expected structural length at time
%     ELw = L_jm/ del_M_A;                                                    % cm, physical length (shell height)

        % time-length-soma weight, Paracas Bay, data from Arturo's thesis
    % ------------------------------    
    %Oxygen limitation  -- hypoxia events
    % oxygen saturation per hour, temperature data per hour
    s_M = L_j/ L_b;                                                         % -, acceleration factor at f    
    tT = temp.tLC3;                                                        % K, temperature from data
    L0 = Lw0.tLC3 * del_M_A;                                                % cm, initial structural length defined in data
    f = f_Paracas;                                                          % -, scaled function response used to compute initial condition and to get the dynamic
    E0 = f * E_m * L0^3;                                                    % J, initial amount of energy in reserves
    VE0  = [L0^3; E0];                                                      % [cm^3; J], initial values of state variables (volume and reserve)
    [t_sort_tLC3, ~ , it_sort_tLC3] = unique(tLC3(:,1),'sorted');           %take the time of the data by ascending order, with unique values
    S_O2c = 40;                                                             % 40% as the oxygen threshold, the organism reggulatory capacity.
    S_O2 = Osat.tLC3;                                                 % %, oxygen saturation in the environment, from data
    
    [~, VE] =  ode45(@dget_VE_hypo, t_sort_tLC3, VE0,[], f, tT, S_O2, par, cPar, s_M, S_O2c);  % resolution of the dynamic with differential equations
    Lw   = (VE(:,1).^(1/3)) / del_M_A;                                      % cm, physical shell height
    ELC3 = Lw(it_sort_tLC3);                                                % reconstruction, physical shell height for data time
    
    %same thing but just use the time of weight data set to be able to
    %custom the plots after
    [t_sort_tWsC3, ~ , it_sort_tWsC3] = unique(tWsC3(:,1),'sorted'); 
    [~, VE] =  ode45(@dget_VE_hypo, t_sort_tWsC3, VE0,[], f, tT, S_O2, par, cPar, s_M, S_O2c);
    
    W_s = d_V * (((VE(:,1).^(1/3))).^3 + VE(:,2) .* w_E / mu_E / d_E);            % g, somatic dry weight
    EWC3 = W_s(it_sort_tWsC3);                                               % reconstruction, somatic dry weight for data time
      
    
%     %time_length larval stage
%     % ------------------------------
%     %keep f = 1, because continuous feeding in the experiment
%     [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
%     t = tLlarv(:,1);                                                        % d, time with fertilisation at zero, from data
%     kT_M = k_M * TC_tLlarv; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;       % 1/d, rates
%     L_bj = L_b * exp(t * rT_j/ 3);                                          % cm, expected structural length at time
%     ELw_larv = L_bj/ del_M_larv_A;                                        % cm, physical length (shell height)


%     % length-weight
%     % ------------------------------
%     %/!\ E_R not considered, here juveniles so before the puberty
%     % -- before metamorphosis
%     EWd_larv = d_V * ((LWlarv(:,1)*del_M_larv_A ).^3  *(1 + f_Chile * w )) ;  % g, dry weight, without reproduction buffer
%     % -- after metamorphosis
%     EWd = d_V * ((LWmeta(:,1)*del_M_A ).^3  *(1 + f_Chile * w )) ;              % g, dry weight, without reproduction buffer

    % pack to output
     prdData.tLC3 = ELC3;
     prdData.tWsC3 = EWC3;
%     prdData.tL = ELw;
%     prdData.tLlarv = ELw_larv;
%     prdData.LWmeta = EWd;
%     prdData.LWlarv = EWd_larv;
    % prdData.LW = EWd;
end

%modified on 20/09/2023 - Eline Le Moan 
%Add the impact of hypoxia, as in Aguirre-Velarde et al., 2019.
function dVE = dget_VE_hypo(t, VE, f, tT, tS, p, c, s_M, S_O2c) 
    % INPUTS: time, initial values (volume and reserve), food, temperature, oxygen saturation, primary parameters, compound parameters, acceleration factor
    % OUTPUTS: derivative of volume and reserve
    % modified from function in Pecten maximus
    
        V = VE(1);     % cm^3, structural volume
        E = VE(2);     % J, reserve
        
        % Temperature
        if length(tT) > 1
            tT_int = interp1(tT(:,1),tT(:,2),t,'linear', 'extrap');  % interpolate the right temperature                
        else
            tT_int = tT; % get the temperature value
        end
        TC = tempcorr(tT_int, p.T_ref, p.T_A); % -, temperature correction factor
        
        % Oxygen saturation
        % after try to put a condition to be able to consider this function
        % even without oxygen saturation (if length(tS) ==1)
        if length(tS) > 1
            tS_int = interp1(tS(:,1),tS(:,2),t,'linear', 'extrap');  % interpolate the right oxygen saturation                
        else
            tS_int = tS; % get the oxygen saturation value
        end
        if tS_int < S_O2c
            C_DO = tS_int / S_O2c;
        else 
            C_DO = 1;
        end

        % Shape correction function applies to surface-area specific assimilation and energy conductance
        % all parameters with time in dimension need to be corrected by the temperature
        p.v    = s_M * p.v * TC;            % cm d^-1, conductance corrected
        c.p_Am = s_M * c.p_Am * TC;         % J cm^-2 d^-1, maximum specific assimilation rate corrected
        p.p_M  = p.p_M * TC;                % J cm^-3 d^-1, volume-specific somatic maintenance rate

        % Rates
        pA  = C_DO * f * c.p_Am * V^(2/3);                                            % J d^-1, assimilation rate
        pC  = C_DO * E * (p.E_G * p.v/ V^(1/3)  + p.p_M)/ (p.kap * E/V + p.E_G);       % J d^-1, mobilisation rate
        pS = p.p_M * V + p.p_T * V^(2/3);                                      % J d^-1, somatic maintenance rate
        pG = p.kap * pC - pS;                                                % J d^-1, growth rate

        % Differential equations of state variables  
        dE = pA - pC;               % J d^-1, reserve dynamic
        dV = pG / p.E_G;           % J d^-1, volume dynamic

        %test to see if V decreases
        if dV < 0
            dV = 0;
        end
        % Pack output 
        dVE = [dV; dE]; 
end

%modified on 19/09/2023 - Eline Le Moan
function dVE = dget_VE(t, VE, f, tT, p, c, s_M) 
    % INPUTS: time, initial values (volume and reserve), food, temperature, primary parameters, compound parameters, acceleration factor
    % OUTPUTS: derivative of volume and reserve
    % modified from function in Pecten maximus
    
        V = VE(1);     % cm^3, structural volume
        E = VE(2);     % J, reserve
        if length(tT) > 1
            tT_int = interp1(tT(:,1),tT(:,2),t,'linear', 'extrap');  % interpolate the right temperature                
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