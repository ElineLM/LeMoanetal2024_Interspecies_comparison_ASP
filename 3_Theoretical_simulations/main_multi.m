%------------------------------------------------------
% Objectives: DEB model simulations for several individuals
% Need to modify init and indiv files for simulations and after run the
% whole main_multi to get results and graphs. 
%
% Outputs:
%         simu: structure with individual features and model predictions
%         graphs of state variables, observable variables
%
% calls:      init.m
%             indiv.m
%             observables.m
%             plots_obs_multi.m
%             plots_SV_multi.m

% 2024/02/06 - Eline Le Moan - based on Laure Pecquerie's codes
% Mascoet Project
%------------------------------------------------------
clear all
close all

global nb_indiv

%% 1 - Initialisation
nb_indiv = 5;       % number of individuals simulated

% in the study, species: Placopecten magellanicus, Nodipecten subnodosus,
% Pecten maximus, Argopecten purpuratus and Mimachlamys varia
colours = ["#E76BF3", "#00BF7D", "#00B0F6", "#F8766D", "#A3A500"];          % colours for graphs
markers = ["v",  "^", "d", "s", "o"];                                       % markers for graphs

%%% the first species is the reference one
ultiHeight = [21, 19, 12, 11, 7];                                           % ultimate length from literature
% enviroT = [8.5, 25, 12, 15, 13]; % Pl, N, P, A, M                         % environmental mean temperature in °C
lifespan = [20*365, 7*365, 20*365, 5*365, 7*365];                           % lifespan from literature

for i = 1:nb_indiv
    allsimu(i).species = init();
    % call parameters from a file, and then here can modify one or several
    % for one or several individuals.
    
    % calculation of ultimate length from parameters chosen
    allsimu(i).species.Licalc = allsimu(i).species.parC.L_m *allsimu(i).species.parStat.s_M / allsimu(i).species.par.del_M;
    % calculation of zoom factor to be used to compare species. From previous Li calculated as reference
    allsimu(i).species.z_species = (ultiHeight(i) * allsimu(i).species.par.del_M) / (ultiHeight(1) * allsimu(1).species.par.del_M) ;
    
    % value to calculate the initial conditions 
    f_init = 1;
%     allsimu(i).species.init.Lp = allsimu(i).species.parStat.L_p * allsimu(i).species.z_species ;
    allsimu(i).species.init.Lp = allsimu(i).species.parStat.L_p  ; % all the same initial length
    %allsimu(i).species.init.Wdp = allsimu(i).species.par.d_V * (allsimu(i).species.init.Lp^3 * (1+allsimu(i).species.parC.w * f_init));
    
    % modification of simulation end time, because consider lifespan of each species
    allsimu(i).species.t_final = lifespan(i);
    % thus modification of simulation time
    allsimu(i).species.time = [allsimu(i).species.t_init : allsimu(i).species.t_final]';
    
    % create the simu structure for calculation for each individual
    simu = allsimu(i).species;
    
    % save the colour and marker type per species 
    allsimu(i).species.colorObs = colours(i);   % colour per species
    allsimu(i).species.marker = markers(i);     % marker type per species    
    
    %% 2 - calcultate state variables
    % call the indiv function to calculate the state variables (E reserve,
    % V structural volume, H maturity, R reproduction buffer)
    % and give a structure with results as vector for time (t) and
    % variables
    simu.tEVHR = indiv(simu);
    allsimu(i).tEVHR = simu.tEVHR;
    
    %% 3 - calculate observable variables
    % from state variable results, call the observables function to
    % calculate observable variables. And save the result in the same
    % structure as before, per species
    simu.obs = observables(simu);
    allsimu(i).obs = simu.obs;
end

%% 4 - plots

%%% Plots of state variables
plots_SV_multi(allsimu)

% %%% Plots of observable variables
plots_obs_multi(allsimu)