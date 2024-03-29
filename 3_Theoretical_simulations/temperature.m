function T = temperature(t, simu)
%------------------------------------------------------
% Objectives: calculate the value of temperature in the environment at each time

% Inputs: 
%             t,      day, time of simulation
%             simu,   structure for species, simulation and enviro information
% Outputs: 
%             T, temperature for simulation in K


% called by:    flux

% 2024/02/06 - Eline Le Moan -  based on Laure Pecquerie's codes
% Mascoet Project
%------------------------------------------------------

% 1 for constant environment, 2 for variable environment
if simu.envT == 1
    T = simu.Tinit;

elseif simu.envT == 2
    T = simu.mean + simu.ampl * sin(2 * pi * (t + simu.shift)/simu.period) ;

end