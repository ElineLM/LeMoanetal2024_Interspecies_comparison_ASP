function simu = init()
%------------------------------------------------------
% Objectives: get all initial values to realise DEB simulation

% Outputs: simu
%
% calls:   one file with parameter values, such as the result file from
% parameter prediction with AmP routine. 


% called by:    main.m

% 2024/02/06 - Eline Le Moan - based on Laure Pecquerie's codes
% Mascoet Project
%------------------------------------------------------

%% Simulation time
simu.t_init = 1;
simu.t_final = 15*365;
simu.time = [simu.t_init : simu.t_final]';


%% Parameters from previous estimations
% dir(fileparts(pwd) + "\" + folderResultsName) % to know which folders and
% files are present in one folder


% get folder name of the results wanted
folderResultsName = "2 - Codes a modifier individuel (5 especes)\Placopecten_magellanicus\results_06-Mar-2024\";
% get file name of the results wanted
fileResultsName = "Newparsinit_dV02_results_Placopecten_magellanicus.mat";

% create the full name of the path to get the right file
resultsPath = fileparts(pwd) + "\" + folderResultsName + fileResultsName;

% open the file
load(resultsPath)

simu.pets = metaData.species;
simu.par = par;
simu.parC = parscomp_st(par);
simu.parStat = statistics_st('abj', simu.par, C2K(20), 1);

% shape coefficient -- could have different name per species
names_fields = fieldnames(simu.par);
index = find(contains(names_fields, 'del_M'));
simu.par.del_M = simu.par.(names_fields{index(1)}); % to get the value of del_M for adult for our species

%% Environmental conditions
% temperature
simu.envT = 1; % constant environment (1) or variable environment (2)

simu.Tinit = C2K(20);
simu.mean = C2K(20);    % averaged temperature in K
simu.ampl = 5;          % difference, no unit
simu.shift = 100;       % shift in time, day
simu.period = 365;      % period, day

% food
simu.envX = 1; % constant environment (1) or variable environment (2)

simu.Xinit = 5;
simu.Xparam = [1460	1.3059
                -0.81167	0.41211
                0.32261	-0.84863
                -0.058888	0.4888
                0.00838	-0.16338
                0.49188	-0.11107
                -0.12596	0.16797
                -0.15009	-0.3274
                0.13013	-0.14936
                0.033949	0.09282
                0.091228	-0.11154];

simu.pressure = 101325;     % Pa (= 1 atm)
simu.ctegaz = 8.3144621;    % J/mol/K

%% Specificity of individuals
%  simu.X_K = 0.5;
 simu.X_K = 0; % to get f = 1 with the Holling response (f = X/(X+X_K))
 simu.t_spawn = 1; % spawn every year - one event


%% Colors for plots
simu.colorSV = "#700404";
simu.colorObs = "#700404";
simu.marker = "o";

end