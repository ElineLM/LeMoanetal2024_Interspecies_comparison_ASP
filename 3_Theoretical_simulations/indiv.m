function tEVHR = indiv(simu)

%------------------------------------------------------
% Objectives: Calculate the initial values and do the integration for state
% variables. In this file we can precise the type of reproduction and all
% other specific characteristic for each individual simulated.
% Here, the same simulation for all species with one spawning per year
% after puberty.

% Inputs:       simu: structure with parameters, and all values needed

% Outputs:      tEVHR: table of time, reserve, structure, maturity and reproduction buffer (results of derivative calculation)
   
% calls:        flux.m

% called by:    main.m

% 2024/02/06 - Eline Le Moan -  based on Laure Pecquerie's codes
% Mascoet Project
%------------------------------------------------------

%% Initialisation
tEVHR = zeros(0,5);
year = 0;

%% Initialisation of state variables
% to start at adult stage, with different values depending on length at puberty (Lp) and EHp
% here only for one individual
f_init = 1;
V_0 = (simu.init.Lp)^3; %L_p given here as structural length in parameters
E_0 = f_init * ((simu.parC.p_Am * simu.z_species) / simu.par.v) * V_0; % no need to correct p_Am and v by temperature and acceleration because ratio between both
EH_0 = simu.par.E_Hp * (simu.z_species^3);
ER_0 = 10^-8; 

simu.EVHR_init = [E_0, V_0, EH_0, ER_0]; % initial values of state variables

%% Call parameters and initial values
EVHR_i = simu.EVHR_init;    % initial values of state variables
t_current = simu.t_init;    % set current time at initial time value

%% Calculation of state variables with spawning events

while t_current < simu.t_final
    year = year + 1;                        % one spawning every year
    t_next = simu.t_spawn + 365 * year;        % for next spawning event, and conversion from year to days
    if t_next > simu.t_final
        t_next = simu.t_final;              % in case the calculation gives higher value
    end
    time = [t_current : t_next]';           % for each year simulated, calculate the state variables
    [t, EVHR] = ode45(@(t, EVHR)flux(time, EVHR, simu), time, EVHR_i);
    tEVHR = [tEVHR; [t, EVHR]];
    E_Hcurrent = tEVHR(end,4);              % get the last value of E_H from simulation
    if E_Hcurrent >= simu.par.E_Hp * simu.z_species^3
        tEVHR(end,5) = 0;                   % check if we are after puberty, if yes, spawning, thus empty E_R = 0
    end
    EVHR_i = tEVHR(end,2:5);
    t_current = tEVHR(end,1);
end

end