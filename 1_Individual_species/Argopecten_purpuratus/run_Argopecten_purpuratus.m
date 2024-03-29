close all; 
global pets 

pets = {'Argopecten_purpuratus'}; 
% global refPets
% refPets = {'Argopecten_purpuratus', 'Nodipecten_subnodosus','Mimachlamys_varia', 'Pecten_maximus','Placopecten_magellanicus'};
check_my_pet(pets); 

estim_options('default'); 
estim_options('max_step_number', 5e2); 
estim_options('max_fun_evals', 5e3); 

estim_options('pars_init_method', 2); 
% estim_options('results_output', 3); 
estim_options('results_output', 3); 
estim_options('method', 'nm'); 

estim_pars; 

load('results_Argopecten_purpuratus.mat')
other_param = statistics_st('abj', par, C2K(18), 1);
other_param.p_Am
%return
% 
%%
estim_options('pars_init_method', 1);
estim_pars;
%
% estim_options('pars_init_method', 1);
% estim_pars;
%
estim_options('pars_init_method', 1);
estim_pars;
%
estim_options('pars_init_method', 1);
estim_options('results_output', 3);
estim_pars;
%
return
%%
% load('results_Argopecten_purpuratus.mat')
% other_param = statistics_st('abj', par, C2K(18), 1);
% other_param.s_M


% Save results to working folder
% commandes modifiées d'après Paulo pour Nodipecten subnodosus
figures = dir('*.png');
html    = dir('*.html');
results = dir('*.mat');

results_file = {figures.name,...
    html.name,...
    results.name,...
    };
timeStamp = char(datetime('today'));            % get the date
saveDir   = ['results_', timeStamp, '/'];       % create the name of the folder with results and the date
trialname = 'fixkappa_dV02_final_';                         % Here change the name of the trial
mkdir(saveDir);

for i = 1:length(results_file)
    copyfile(results_file{i},[saveDir, trialname, results_file{i}])
end