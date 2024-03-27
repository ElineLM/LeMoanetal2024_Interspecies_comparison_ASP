function [prdData, info] = predict_Pecten_maximus(par, data, auxData)

    % unpack par, data, auxData
   %%% calculation for body size scaling relationship
    par.E_Hp = par.E_Hp_ref*((par.z_species)^3);
    cPar = parscomp_st(par);
    cPar.p_Am = (par.z_ref * par.z_species) * par.p_M/ par.kap;   % J/d.cm^2, {p_Am} spec assimilation flux; the expression for p_Am is multiplied also by L_m^ref = 1 cm, for units to match. 
    vars_pull(par);   vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

    % del_M max fixed at 1.3 as the shape coeff for a sphere
    filterChecks = ...% f contrained to not be larger than 1.2 or negative
        f_EV > 1.2 || f_EV < 0 || ...; %
        f_B > 1.2  || f_B < 0 ;%|| ...; %
%         del_M_N < del_M_P || del_M_P < del_M_A || del_M_A < del_M_Pl || del_M_Pl < del_M_M;
                   
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
    TC_EV = tempcorr(temp.tLl_EV, T_ref, T_A);
    %   TC_tLlarvae = tempcorr(temp.tLlarvae, T_ref, T_A);
  
    % -----------------------------------------------------------------------
    % zero-variate data
    % -----------------------------------------------------------------------
    % life cycle
    pars_tj = [g k l_T v_Hb v_Hj v_Hp];
    [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

    pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
    U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
    E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap and v constant

    % birth
    L_b = L_m * l_b;                          % cm, structural length at birth
    Lw_b = L_b/ del_M_larv_P;                 % cm, shell height at birth
    Wd_b = L_b^3 * d_V * (1 + f * w);         % g, dry weight at birth
    aT_b = t_b/ k_M/ TC_ab;                   % d, age at birth 

    % metam
    L_j = L_m * l_j;                          % cm, structural length at metam
    Lw_j = L_j/ del_M_P;                      % cm, shell height at metam
    Wd_j = L_j^3 * d_V * (1 + f * w);         % g, dry weight at metam
    tT_j = (t_j - t_b)/ k_M/ TC_tj;           % d, time since birth at metam

    % puberty 
    L_p = L_m * l_p;                          % cm, structural length at puberty 
    Lw_p = L_p/ del_M_P;                       % cm, shell height at puberty 
    Wd_p = L_p^3 * d_V * (1 + f * w);         % g, dry weight at puberty 
    tT_p = (t_p - t_b)/ k_M/ TC_tp;           % d, time since birth at puberty 

    % ultimate
    L_i = L_m * l_i;                          % cm, ultimate structural length
    Lw_i = L_i/ del_M_P;                        % cm, ultimate shell height
    Wd_i = L_i^3 * d_V * (1 + f * w);         % g, ultimate dry weight (without reproduction buffer)

    % reproduction
%     pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector at T_ref
%     RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % #/d, ultimate reproduction rate at T
    pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp];   % compose parameter vector at T
    LR = L_R.R_L;                                                   % length at which reproductive effort is given
    RT_L = TC_RL * reprod_rate_j(LR, f, pars_R);                    % #/d, ultimate reproduction rate at T

    % life span
    pars_tm = [g; l_T; h_a/ k_M^2; s_G];      % compose parameter vector at T_ref
    t_m = get_tm_s(pars_tm, f, l_b);          % -, scaled mean life span at T_ref
    aT_m = t_m/ k_M/ TC_am;                   % d, mean life span at T

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

    % -----------------------------------------------------------------------
    % uni-variate data
    % -----------------------------------------------------------------------

%     % time-length (celtic sea)
%     % ------------------------------
%     s_M = L_j/ L_b;                   % -, acceleration factor at f
%     time = tL1(:,1);                  % d, from data
%     tT = temp.tL1;                    % K, temperature from data
%     L0 = Lw0.tL1 * del_M_P;           % cm, initial structural length
%     f0 = f_tL1;                       % -, initial scaled function response is used to compute initial condition of the individual
%     E0 = f0 * E_m * L0^3;             % J, initial amount of energy in reserves
%     VE0 = [L0^3; E0];                 % [cm^3, J], initial values of the state variables (structure and reserve)
%     [~, VE] =  ode45(@dget_VE, time, VE0,[],f, tT, par, cPar, s_M);   % 4-1 vector [J, cm, J, J] with initial conditions [L0; E0; EH0; ER0]
%     L   = VE(:,1).^(1/3);             % cm, structural length
%     ELw1 = L/del_M_P;                 % cm, physical shell height


    % time-length, bay of brest, France
    % ------------------------------
    s_M = L_j/ L_b;                 % -, acceleration factor at f
    time = tLl_B(:,1);                % d, from data
    tT = temp.tLl_B;                  % K, temperature from data
    L0 = Lw0.tLl_B * del_M_P;         % cm, initial structural length
    f = f_B;                      % -, initial scaled function response is used to compute initial condition of the individual
    E0 = f * E_m * L0^3;            % J, initial amount of energy in reserves
    VE0 = [L0^3; E0];               % [cm^3, J], initial values of the state variables (structure and reserve)
    [~, VE] =  ode45(@dget_VE, time, VE0,[],f, tT, par, cPar, s_M);
    L   = VE(:,1).^(1/3);           % cm, structural length
    ELw_B = L/del_M_P;                 % cm, physical shell height

%     % % time-length, traena, Norway
%     % ------------------------------ 
%     s_M = L_j/ L_b;                   % -, acceleration factor at f
%     time = tL3(:,1);                  % d, from data
%     tT = temp.tL3;                    % K, temperature from data
%     L0 = Lw0.tL3 * del_M_P;           % cm, initial structural length
%     f = f_tL3;                        % -, initial scaled function response is used to compute initial condition of the individual
%     E0 = f * E_m * L0^3;              % J, initial amount of energy in reserves
%     VE0 = [L0^3; E0];                 % [cm^3, J], initial values of the state variables (structure and reserve)
%     [~, VE] =  ode45(@dget_VE, time, VE0,[],f, tT, par, cPar, s_M); 
%     L   = VE(:,1).^(1/3);             % cm, structural length
%     ELw3 = L/del_M_P;                 % cm, physical shell height

    % time-length, Bay of Brest, EVECOS
    % ------------------------------
    t = tLl_EV(:,1);                                                            % d, time since birth
    f = f_EV;                                                               % -, mean scaled func response for tL
    [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);
    L_b = l_b * L_m; L_j = l_j * L_m; L_i = l_i * L_m;                      % cm, struct. length at birth, at metamorphosis and ultimate
    kT_M = k_M * TC_EV; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M;           % 1/d, rates
    tT_j = (t_j - t_b)/ kT_M;                                               % only after metamorphosis here in the data  % d, time end of acceleration
    L_jm = L_i - (L_i - L_j) * exp( - rT_B * (t(t >= tT_j) - tT_j));        % cm, structural length after metamorphosis
    ELw_EV = L_jm/ del_M_P;                                                      % cm, physical length (shell height)
    
%     % time-dry weight, Bay of Brest, EVECO (retrieved from graph in
%     % Lavaud's thesis) (18/09/2023)
%     % ------------------------------
%     s_M = L_j/ L_b;                 % -, acceleration factor at f    
%     time = tWd(:,1);                % d, from data
%     tT = temp.tWd;                  % K, temperature from data
%     L0 = Lw0.tWd * del_M_P;         % cm, initial structural length
%     f = f_tL;                       % -, initial scaled function response is used to compute initial condition of the individual
%     E0 = f * E_m * L0^3;            % J, initial amount of energy in reserves
%     VE0  = [L0^3; E0];              % [cm^3, J], initial values of the state variables (structure and reserve)
%     [~, VE] =  ode45(@dget_VE, time, VE0, [], f, tT, par, cPar, s_M); 
%     Lw   = (VE(:,1)^(1/3))/del_M_P;                                   % cm, physical shell height
%     EtWd = d_V * (VE(:,1).^3 + VE(:,2) .* w_E / mu_E / d_E);         % g, somatic dry weight
     
%     % Length-weight, Galicia
%     % ------------------------------
%     EWd_Gal = d_V * ((LW_Galicia(:,1)*del_M_P ).^3  *(1 + f_Gal * w )) ;      % g, dry somatic weight
%     % consider the shell height and no reproduction buffer in the weight

%     % time-length for early stages 
%     % ------------------------------
%     % Data from the Tinduff hatchery, controled temperature and food
%     [t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_tj(pars_tj, f_Tinduff); 
%     %f_Tinduff, the functional response for all data sets coming from Tinduff hatchery
%     %%% Larvae 
%     rT_j = TC_tLlarvae * rho_j * k_M;                 % exponential growth rate, corrected with the temperature
%     time_bj = tLlarvae(:,1);                          % days, time between birth and metamorphosis
%     Lw_b = l_b * L_m/ del_M_larv_P;                   % cm, physical length at birth
%     ELw_bj = Lw_b * exp(time_bj * rT_j/3);            % cm, physical length from birth to metamrorphosis (shell length)


    % pack to output
    prdData.tLl_EV = ELw_EV;
%     prdData.tL1 = ELw1;
    prdData.tLl_B = ELw_B;
%     prdData.tWd = EtWd;
%     prdData.tL3 = ELw3;
%     prdData.LW_Galicia = EWd_Gal; 
%     prdData.tLlarvae = ELw_bj;
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