% close all
global pets

% species names for the group (alphabetical order)
pets = {'Argopecten_purpuratus', 'Mimachlamys_varia', 'Nodipecten_subnodosus', 'Pecten_maximus', 'Placopecten_magellanicus'};

estim_options('default');             % use Nead-Melder simplex method and filter

estim_options('max_step_number',5e2); 
% estim_options('max_step_number',10e2); 
estim_options('max_fun_evals',5e3);   

estim_options('pars_init_method', 2);
estim_options('results_output',3);
estim_options('method', 'nm');

estim_pars; 


%% continue the estimation
% 
estim_options('pars_init_method', 1);
estim_pars; 

return 
%(2 continuous estimations)
estim_options('pars_init_method', 1);
estim_pars; 

estim_options('pars_init_method', 1);
estim_pars; 

estim_options('pars_init_method', 1);
estim_pars; 

estim_options('max_step_number',10e2); 
estim_options('pars_init_method', 1);
estim_pars; 

estim_options('pars_init_method', 1);
estim_options('results_output', 3);
estim_pars; 

% return

% %%
% load('results_group.mat')
% other_param = statistics_st('abj', par, C2K(18), 1);
% other_param.p_Am

%% Save results to working folder
% modified lines from Paulo Lagos for Nodipecten subnodosus
figures = dir('*.png');
html    = dir('*.html');
results = dir('*.mat');

results_file = {figures.name,...
    html.name,...
    results.name,...
    };
timeStamp = char(datetime('today'));            % get the date
saveDir   = ['results_', timeStamp, '/'];       % create the name of the folder with results and the date
trialname = 'estim_f1init_run3_';                         % Here change the name of the trial
mkdir(saveDir);


for i = 1:length(results_file)
    copyfile(results_file{i},[saveDir, trialname, results_file{i}])
end


return 

% %% Save parameter values into similar file for all species
% % in csv files
% trialnameSave = 'kappa_bivalves_fixed_2_';
% load([trialnameSave, 'results_group.mat'])
% temp_par = rmfield(par, 'free');
% table_par = struct2table(temp_par);
% writetable(table_par, [timeStamp,'_', trialnameSave, 'parameters.csv'])

