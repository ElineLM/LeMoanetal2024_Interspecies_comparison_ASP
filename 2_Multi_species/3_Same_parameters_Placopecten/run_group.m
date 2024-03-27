% close all
global pets

% species names for the group (alphabetical order)
pets = {'Argopecten_purpuratus', 'Mimachlamys_varia', 'Nodipecten_subnodosus', 'Pecten_maximus', 'Placopecten_magellanicus'};

estim_options('default');             % use Nead-Melder simplex method and filter

estim_options('max_step_number',5e2);
estim_options('max_fun_evals',5e3);   

estim_options('pars_init_method', 2);
estim_options('results_output',3);
estim_options('method', 'no');

estim_pars; 

return

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
trialname = 'trial_';  % Here change the name of the trial
mkdir(saveDir);


for i = 1:length(results_file)
    copyfile(results_file{i},[saveDir, trialname, results_file{i}])
end