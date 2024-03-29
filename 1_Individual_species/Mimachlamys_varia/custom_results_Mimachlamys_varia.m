%% custom_results_my_pet
% presents results of univariate data graphically in a customized way

%%
function custom_results_Mimachlamys_varia(par, metaPar, data, txtData, auxData)
%par, metaPar, txtPar, data, auxData, metaData, txtData, weights
  % created by Starrlight Augustine, Dina Lika, Bas Kooijman, Goncalo Marques and Laure Pecquerie 2015/04/12
  % modified 2015/08/25
  % modified 2018/09/06 by Nina Marn
  
  %% Syntax
  % <../custom_results_template.m *custom_results_my_pet*>(par, metaPar, txtData, data, auxData)
  
  %% Description
  % present customized results of univariate data
  %
  % Inputs:
  %
  % * par: structure with parameters (see below)
  % * metaPar: structure with field T_ref for reference temperature
  % * txt_data: text vector for the presentation of results
  % * data: structure with data
  % * auxData: structure with temperature data and potential food data
  
  %% Remarks
  % * A template named 'custom_results_template' is available in 'pet' folder of DEBtool_M:  
  % Replace '_template' in the function name with 'my_pet' to use with my_pet templates 
  % * Modify to select and plot uni-variate data for your entry: copy to  folder of your species, 
  %     replacing 'template' (or 'my_pet') with the Latin name of your species,
  %     and template data with the entry-specific data you wish to plot
  % * Once named appropriately, this function will be called automatically by 
  %     <results_pets.html *results_pets*> function of DEBtool_M when running the <run_my_pet.html *run*> file
  
  try
      
  % get predictions
  data2plot = data;              % copy data to Prd_data
%   t = linspace(150,1200, 262)';   % set independent variable
%   L = linspace(1, 6, 262)'; % set independent variable
  t = linspace(min(data2plot.tL19SA(:,1)), max(data2plot.tL19SA(:,1)), 5000)';   % set independent variable
  L = linspace(min(data2plot.LW19SA(:,1)), max(data2plot.LW19SA(:,1)), 5000)';
  data2plot.tL19SA = t; % overwrite independent variable in tL
  data2plot.tW19SA = t; % overwrite independent variable in tW
  data2plot.LW19SA = L; % overwrite independent variable in LW
  [prdData, info] = predict_Mimachlamys_varia(par, data2plot, auxData);
  
  [stat, txt_stat]  = feval('statistics_st', metaPar.model, par, C2K(20), par.f);
 
  if strcmp(metaPar.model, 'abj')
    fprintf(['\n acceleration factor s_M is ', num2str(stat.s_M), ' \n'])
  end

  % unpack data & predictions
  tL     = data.tL19SA;     % data points first set
  tW     = data.tW19SA;     % data points second set
  LW     = data.LW19SA;     % data points thirt set
  
  EL     = prdData.tL19SA; % predictions (dependent variable) first set
  EW     = prdData.tW19SA; % predictions (dependent variable) second set
  ELW    = prdData.LW19SA; % predictions (dependent variable) third set

  close all % remove existing figures, else you get more and more if you retry

  figure('Name', 'data1') % figure to show results of uni-variate data
  %characteristics for both figures
  set(gca,'FontSize',22, 'Box', 'on')
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperUnits','points'); 
  set(gcf, 'Position',  [120, 120, 500, 400])
  plot(t, EL, 'Color', '#1886ad', 'linewidth', 2)
  hold on
  plot(tL(:,1), tL(:,2), '.', "MarkerEdgeColor", "#940f06", 'MarkerFaceColor', "#940f06", 'markersize', 20)
  xlabel([txtData.label.tL19SA{1}, ', ', txtData.units.tL19SA{1}])%, 'FontSize', 20)
  ylabel([txtData.label.tL19SA{2}, ', ', txtData.units.tL19SA{2}])%, 'FontSize', 20)
  set(findall(gcf,'-property','FontSize'),'FontSize',24)
  print -dpng results_Mimachlamys_varia_01.png
   
   
  figure('Name', 'data2') % figure to show results of uni-variate data
  plot(t, EW, 'Color', '#1886ad', 'linewidth', 2)
  hold on
  plot(tW(:,1), tW(:,2), '.', "MarkerEdgeColor", "#940f06", 'MarkerFaceColor', "#940f06", 'markersize', 20)
  xlabel([txtData.label.tW19SA{1}, ', ', txtData.units.tW19SA{1}])
  ylabel({[txtData.label.tW19SA{2}, ',']; txtData.units.tW19SA{2}})
  set(findall(gcf,'-property','FontSize'),'FontSize',24)
  print -dpng results_Mimachlamys_varia_02.png
  
    figure('Name', 'data3') % figure to show results of uni-variate data
  plot(L, ELW, 'Color', '#1886ad', 'linewidth', 2)
  hold on
  plot(LW(:,1), LW(:,2), '.', "MarkerEdgeColor", "#940f06", 'MarkerFaceColor', "#940f06", 'markersize', 20)
  xlabel([txtData.label.LW19SA{1}, ', ', txtData.units.LW19SA{1}])
  ylabel({[txtData.label.LW19SA{2}, ',']; txtData.units.LW19SA{2}})
  set(findall(gcf,'-property','FontSize'),'FontSize',24)
  print -dpng results_Mimachlamys_varia_03.png
  
  catch
    fprintf('Warning from custom_results_template: this template is meant to replace the default way of presenting results\n');
    fprintf('This file requires case-specific editing, which is not yet done properly\n');
    fprintf('Use the default way of presenting results by removing this file from your current local directory\n');
    return
 end
