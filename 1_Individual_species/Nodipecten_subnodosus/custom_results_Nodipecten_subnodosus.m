%% custom_results_my_pet
% presents results of univariate data graphically in a customized way

%%
function custom_results_Nodipecten_subnodosus(par, metaPar, data, txtData, auxData)
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
  
  
      global pets
   
      n_fields = length(fieldnames(data));
      names_fields = fieldnames(data);
      for j = 1:n_fields
          % length at age
          if contains(names_fields{j}, 'tL') == 1        
              figure
              data2plot = data;
              tLi = data.(names_fields{j});
              t = linspace(min(tLi(:,1)), max(tLi(:,1)), 5000)';
              data2plot.(names_fields{j}) = t;

%               eval(['prdDatai = predict_', pets, '(par, data2plot, auxData);']);
              [prdData, info] = predict_Nodipecten_subnodosus(par, data2plot, auxData);
              tLp = prdData.(names_fields{j});
              tp = data2plot.(names_fields{j});

              plot(tLi(:,1), tLi(:,2), '.', 'color', "#940f06", 'Markersize', 15);
              hold on
              plot(tp, tLp, '--', 'color', '#1886ad', 'linewidth', 2);
              xlabel([txtData.label.(names_fields{j}){1}, ', ', txtData.units.(names_fields{j}){1}]);
              ylabel([txtData.label.(names_fields{j}){2}, ', ', txtData.units.(names_fields{j}){2}])
%               xlabel("Time (day)");
%               ylabel("Shell height (cm)"); 
              set(gca,'Fontsize',18, 'Box', 'on')
              print(['results_Nodipecten_subnodosus', '_', names_fields{j}], '-dpng')
                         
          % weight at age
          elseif contains(names_fields{j}, 'tW') == 1
              figure
              data2plot = data;            
              tWi = data.(names_fields{j});
              t = linspace(min(tWi(:,1)), max(tWi(:,1)), 5000)';
              data2plot.(names_fields{j}) = t;
              
%               eval(['prdDatai = predict_', pets, '(par, data2plot, auxData;']);
              [prdData, info] = predict_Nodipecten_subnodosus(par, data2plot, auxData);

              tWp = prdData.(names_fields{j});
              tp = data2plot.(names_fields{j});
              
              plot(tWi(:,1), tWi(:,2), '.', 'color', "#940f06", 'Markersize', 15);
              hold on
              plot(tp, tWp, '--', 'color', '#1886ad', 'linewidth', 2);
              xlabel([txtData.label.(names_fields{j}){1}, ', ', txtData.units.(names_fields{j}){1}]);
              ylabel([txtData.label.(names_fields{j}){2}, ', ', txtData.units.(names_fields{j}){2}]); 
%               xlabel("Age (day)");
%               ylabel("Tissue weight (g)"); 
              set(gca,'Fontsize',18, 'Box', 'on')
              print(['results_Nodipecten_subnodosus', '_', names_fields{j}], '-dpng')
              
          % weight at length   
          elseif contains(names_fields{j}, 'LW') == 1
              figure
              data2plot = data;
              LWi = data.(names_fields{j});
              Lg = linspace(min(LWi(:,1)), max(LWi(:,1)), 5000)';
              data2plot.(names_fields{j}) = Lg;
%               eval(['prdDatai = predict_', pets, '(par, data2plot, auxData;']);
              [prdData, info] = predict_Nodipecten_subnodosus(par, data2plot, auxData);
              LWp = prdData.(names_fields{j});
              Lg = data2plot.(names_fields{j});
              
              plot(LWi(:,1), LWi(:,2), '.', 'color', "#940f06", 'Markersize', 15);
              hold on
              plot(Lg, LWp, '--', 'color','#1886ad', 'linewidth', 2);
              xlabel([txtData.label.(names_fields{j}){1}, ', ', txtData.units.(names_fields{j}){1}]);
              ylabel([txtData.label.(names_fields{j}){2}, ', ', txtData.units.(names_fields{j}){2}]); 
%               xlabel("Shell height (cm)");
%               ylabel("Tissue weight (g)"); 
              set(gca,'Fontsize',18, 'Box', 'on')
              print(['results_Nodipecten_subnodosus', '_', names_fields{j}], '-dpng')
          end   
          
%           print -dpng results_group.png
      end 
      

    
      
      
      
      
      
      
      
      
%       
%       
%       
%   % get predictions
%   data2plot = data;              % copy data to Prd_data
%   t_RC = linspace(min(data2plot.tLh_RC(:,1)), max(data2plot.tLh_RC(:,1)), 5000)';
%   t_RA = linspace(min(data2plot.tLh_RA(:,1)), max(data2plot.tLh_RA(:,1)), 5000)';   % set independent variable
%   t_AM = linspace(min(data2plot.tLh_AM(:,1)), max(data2plot.tLh_AM(:,1)), 5000)';
%   L_CL = linspace(min(data2plot.CL_LWd(:,1)), max(data2plot.CL_LWd(:,1)), 5000)';
%   t_MA = linspace(min(data2plot.tLh_MA(:,1)), max(data2plot.tLh_MA(:,1)), 5000)';
%   data2plot.tLh_RC = t_RC; % overwrite independent variable in tL Racotta
%   data2plot.tLh_RA = t_RA; % overwrite independent variable in tL Ramirez-Arce
%   data2plot.tWw_RA = t_RA; % overwrite independent variable in tWw Ramirez-Arce
%   data2plot.tLh_AM = t_AM; % overwrite independent variable in tL Arellano-Martinez
%   data2plot.CL_LWd = L_CL; % overwrtie independent variable in LWd Cleareboudt
%   data2plot.tLl_MA = t_MA; % overwrtie independent variable in tLl Maldonado-Amparo
%   [prdData, info] = predict_Nodipecten_subnodosus(par, data2plot, auxData);
%   
%   [stat, txt_stat]  = feval('statistics_st', metaPar.model, par, C2K(20), par.f);
%  
%   if strcmp(metaPar.model, 'abj')
%     fprintf(['\n acceleration factor s_M is ', num2str(stat.s_M), ' \n'])
%   end
% 
%   % unpack data & predictions
%   tLRC      = data.tLh_RC;
%   tLRA      = data.tLh_RA;     % data points first set
%   tWRA      = data.tWw_RA;     % data points second set
%   tLAM      = data.tLh_AM;     % data points third set
%   LWCL      = data.CL_LWd;     % data points fourth set
%   tLMA      = data.tLl_MA;
%   
%   ELRC      = prdData.tLh_RC
%   ELRA      = prdData.tLh_RA; % predictions (dependent variable) first set
%   EWRA      = prdData.tWw_RA; % predictions (dependent variable) second set
%   ELAM      = prdData.tLh_AM; % predictions (dependent variable) third set
%   EWCL      = prdData.CL_LWd; % predictions (dependent variable) fourth set
%   ELMA      = prdData.tLl_MA;
% 
%   close all % remove existing figures, else you get more and more if you retry
% 
%   figure('Name', 'data_RA_tL') % figure to show results of uni-variate data
%   %characteristics for both figures
%   set(gca,'FontSize',22, 'Box', 'on')
%   set(gcf,'PaperPositionMode','manual');
%   set(gcf,'PaperUnits','points'); 
%   set(gcf, 'Position',  [5, 5, 700, 600])
%   plot(t_RA, ELRA, 'Color', '#1886ad', 'linewidth', 2)
%   hold on
%   plot(tLRA(:,1), tLRA(:,2), '.', "MarkerEdgeColor", "#940f06", 'MarkerFaceColor', "#940f06", 'markersize', 20)
%   xlabel([txtData.label.tLh_RA{1}, ', ', txtData.units.tLh_RA{1}])%, 'FontSize', 20)
%   ylabel([txtData.label.tLh_RA{2}, ', ', txtData.units.tLh_RA{2}])%, 'FontSize', 20)
%   set(findall(gcf,'-property','FontSize'),'FontSize',24)
%   print -dpng results_Nodipecten_subnodosus_01.png
%    
%    
%   figure('Name', 'data_RA_tWw') % figure to show results of uni-variate data
%   plot(t_RA, EWRA, 'Color', '#1886ad', 'linewidth', 2)
%   hold on
%   set(gcf,'PaperPositionMode','manual');
%   set(gcf,'PaperUnits','points'); 
%   set(gcf, 'Position',  [250, 200, 700, 600])
%   plot(tWRA(:,1), tWRA(:,2), '.', "MarkerEdgeColor", "#940f06", 'MarkerFaceColor', "#940f06", 'markersize', 20)
%   xlabel([txtData.label.tWw_RA{1}, ', ', txtData.units.tWw_RA{1}])
%   ylabel([txtData.label.tWw_RA{2}, ', ', txtData.units.tWw_RA{2}])
%   set(findall(gcf,'-property','FontSize'),'FontSize',24)
%   print -dpng results_Nodipecten_subnodosus_02.png
%   
%   figure('Name', 'data_AM_tL') % figure to show results of uni-variate data
%   %characteristics for both figures
%   plot(tAM, ELAM, 'Color', '#1886ad', 'linewidth', 2)
%   hold on
%   set(gcf,'PaperPositionMode','manual');
%   set(gcf,'PaperUnits','points'); 
%   set(gcf, 'Position',  [900, 200, 700, 600])
%   plot(tLAM(:,1), tLAM(:,2), '.', "MarkerEdgeColor", "#940f06", 'MarkerFaceColor', "#940f06", 'markersize', 20)
%   xlabel([txtData.label.AM_tLh{1}, ', ', txtData.units.AM_tLh{1}])%, 'FontSize', 20)
%   ylabel([txtData.label.AM_tLh{2}, ', ', txtData.units.AM_tLh{2}])%, 'FontSize', 20)
%   set(findall(gcf,'-property','FontSize'),'FontSize',24)
%   print -dpng results_Nodipecten_subnodosus_03.png
%    
%     
%   figure('Name', 'data_CL_LWd') % figure to show results of uni-variate data
%   %characteristics for both figures
%   plot(LCL, EWCL, 'Color', '#1886ad', 'linewidth', 2)
%   hold on
%     set(gcf,'PaperPositionMode','manual');
%   set(gcf,'PaperUnits','points'); 
%   set(gcf, 'Position',  [1000, 5, 700, 600])
%   plot(LWCL(:,1), LWCL(:,2), '.', "MarkerEdgeColor", "#940f06", 'MarkerFaceColor', "#940f06", 'markersize', 20)
%   xlabel([txtData.label.CL_LWd{1}, ', ', txtData.units.CL_LWd{1}])%, 'FontSize', 20)
%   ylabel([txtData.label.CL_LWd{2}, ', ', txtData.units.CL_LWd{2}])%, 'FontSize', 20)
%   set(findall(gcf,'-property','FontSize'),'FontSize',24)
%   print -dpng results_Nodipecten_subnodosus_04.png
%     
%   catch
%     fprintf('Warning from custom_results_template: this template is meant to replace the default way of presenting results\n');
%     fprintf('This file requires case-specific editing, which is not yet done properly\n');
%     fprintf('Use the default way of presenting results by removing this file from your current local directory\n');
%     return
%  end
