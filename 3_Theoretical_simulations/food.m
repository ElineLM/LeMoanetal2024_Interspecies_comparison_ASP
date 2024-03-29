function X = food(t, simu)
%------------------------------------------------------
% Objectives: calculate the value of food in the environment at each time


% Inputs: 
%             t,      day, time of simulation
%             simu,   structure for species, simulation and enviro information
% Outputs: 
%             X,    food for simulation 


% called by:    flux

% 2024/02/06 - Eline Le Moan - based on codes from Laure Pecquerie
% Mascoet Project
%------------------------------------------------------

% 1 for constant environment, 2 for variable environment
if simu.envX == 1
    X = simu.Xinit;

elseif simu.envX == 2
    X = fnfourier(t, simu.Xparam);

end