close all; 
global pets 

pets = {'Mimachlamys_varia'}; 
check_my_pet(pets); 

estim_options('default'); 
estim_options('max_step_number', 5e2);
estim_options('max_fun_evals', 5e3);

estim_options('pars_init_method', 2); 
estim_options('results_output', 3);  
estim_options('method', 'nm'); 



estim_pars; 
load('results_Mimachlamys_varia.mat');
other_param = statistics_st('abj', par, C2K(18), 1); %to have other parameters like s_M
% disp('s_M = ')
% disp(other_param.s_M)
% disp('p_Am = ')
disp(other_param.p_Am)
% disp('E_0 = ')
% disp(other_param.E_0) % current estimation E_0 = 3.10^-3 J cf Eline


% return

%% Continue the estimation
estim_options('pars_init_method', 1); 
% estim_options('results_output', 3);  
estim_pars;
return
%% Save results to working folder
% commandes modifiées d'après Paulo Lagos pour Nodipecten subnodosus
figures = dir('*.png');
html    = dir('*.html');
results = dir('*.mat');

results_file = {figures.name,...
    html.name,...
    results.name,...
    };
timeStamp = char(datetime('today'));            % get the date
saveDir   = ['results_', timeStamp, '/'];       % create the name of the folder with results and the date
trialname = 'dV02_final_';                         % Here change the name of the trial
mkdir(saveDir);

for i = 1:length(results_file)
    copyfile(results_file{i},[saveDir, trialname, results_file{i}])
end