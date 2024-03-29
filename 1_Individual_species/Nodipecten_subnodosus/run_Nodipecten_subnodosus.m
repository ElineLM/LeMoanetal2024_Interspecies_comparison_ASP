close all; 
global pets 

pets = {'Nodipecten_subnodosus'}; 
global refPets
refPets = {'Pecten_maximus','Argopecten_purpuratus','Mimachlamys_varia', 'Placopecten_magellanicus', 'Nodipecten_subnodosus'};
check_my_pet(pets); 


estim_options('default'); 
estim_options('max_step_number',5e2); 
estim_options('max_fun_evals',5e3); 


estim_options('pars_init_method', 2);
estim_options('results_output', 3);
estim_options('method', 'nm'); 
estim_pars;
return
%%
estim_options('pars_init_method', 1);
estim_pars; 

estim_options('pars_init_method', 1);
estim_pars; 
% 
% estim_options('pars_init_method', 1);
% estim_pars; 

estim_options('pars_init_method', 1);
estim_options('results_output', 3);
estim_pars;

load('results_Nodipecten_subnodosus.mat');
other_param = statistics_st('abj', par, C2K(20), 1); %to have other parameters like s_M
disp('s_M = ')
disp(other_param.s_M)
disp('p_Am = ')
disp(other_param.p_Am)
disp('E_0 = ')
disp(other_param.E_0)

return 
%% continue the estimation 
% estim_options('pars_init_method', 1);
% estim_pars; 
%
% estim_options('pars_init_method', 1);
% estim_options('results_output', 3);
% estim_pars; 
%
% other_param = statistics_st('abj', par, C2K(20), 1); %to have other parameters like s_M
% disp('s_M = ')
% disp(other_param.s_M)
% disp('p_Am = ')
% disp(other_param.p_Am)
% disp('E_0 = ')
% disp(other_param.E_0) 

%% Save results to working folder
% script modified from Paulo Lagos
figures = dir('*.png');
html    = dir('*.html');
results = dir('*.mat');

results_file = {figures.name,...
    html.name,...
    results.name,...
    };
timeStamp = char(datetime('today'));            % get the date
saveDir   = ['results_', timeStamp, '/'];       % create the name of the folder with results and the date
trialname = 'dV02_fixTA_w2Lp_6runs_';                         % Here change the name of the trial
mkdir(saveDir);

for i = 1:length(results_file)
    copyfile(results_file{i},[saveDir, trialname, results_file{i}])
end

