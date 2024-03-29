%% custom_results_group
% presents results of univariate data graphically in a customized way for a multi-species case

%%
function custom_results_group(parGrp, metaPar, data, txtData, auxData)
% created by Bas Kooijman 2019/03/11
  
  %% Syntax
  % <../custom_results_group.m *custom_results_group*>(parGrp, metaPar, txtData, data, auxData)
  
  %% Description
  % presents customized results of univariate data in multispecies-plots
  %
  % Inputs:
  %
  % * parGrp: structure with parameters ar specified in the pars_init_grp.m file
  % * metaPar: structure with field T_ref for reference temperature
  % * txt_data: text vector for the presentation of results
  % * data: structure with data
  % * auxData: structure with temperature data and potential food data
  
  %% Remarks
  % if a file with this name is present in the directory of the group, it will suppress other graphical output
  
  try 
  % convert all inputs to structures, where the first field names are the pets (ignorging auxData) 
  parGrp = rmfield(parGrp,'free'); par = parGrp2Pets(parGrp); pets = fieldnames(par); n_pets = length(pets);
  %model = metaPar.model; cov_rules = metaPar.cov_rules; metaPar = rmfield(metaPar, {'model','cov_rules'});

  parFlds = fieldnames(parGrp); n_parFlds = length(parFlds); shareTxt = 'share: '; diffTxt = 'diff: ';
  for i = 1:n_parFlds
    if strcmp(parFlds{i}, 'd_X')
      break
    end
    if length(parGrp.(parFlds{i})) == 1
      shareTxt = [shareTxt, parFlds{i}, ', ']; 
    else
      diffTxt = [diffTxt, parFlds{i}, ', '];
    end
  end
  shareTxt(end - 0:1) = ''; diffTxt(end - 0:1) = '';
   
   
  close all % remove existing figures, else you get more and more if you retry
%   figure % figure to show results of uni-variate data
%   hold on % prepare for multiple plots in one figure

  
%   t = linspace(0, 6000, 15000)';   % set independent variable for predictions (knowing that this applies to all tL-plots)
%   colors = {[1, 0, 0], [1, 0, 1], [0, 0, 1], [0, 0, 0], [1, 1, 0]}; % colors for pet 1,.,4 and 5
  colors = {'#a10505', '#e38f19', '#1b5780','#22871c', '#ab0777'}; % Argopecten, Mimachlamys, Nodipecten, Pecten, Placopecten
%   colors = { '#a10505', '#e38f19','#ab0777'}; %Argopecten, Mimachlamys, Placopecten

  for i = 1:n_pets
      n_fields = length(fieldnames(data.(pets{i})));
      names_fields = fieldnames(data.(pets{i}));
      for j = 1:n_fields
          % length at age
          if contains(names_fields{j}, 'tL') == 1
              if sum(contains(names_fields, 'tL')) > 1
                  names_fields_tL = names_fields(contains(names_fields, 'tL'));
                  figure
                  data2plot = data.(pets{i});
                  for k = 1:length(names_fields_tL)               
                      tLi = data.(pets{i}).(names_fields_tL{k});
                      t = linspace(min(tLi(:,1)), max(tLi(:,1)), 5000)';
                      data2plot.(names_fields_tL{k}) = t;

        %               eval(['prdDatai = predict_', pets{i}, '(par.(pets{i}), data.(pets{i}), auxData.(pets{i}));']);
                      eval(['prdDatai = predict_', pets{i}, '(par.(pets{i}), data2plot, auxData.(pets{i}));']);
        %               tLp = unique(prdDatai.(names_fields{j}));
        %               tp = linspace(min(tLi(:,1)), max(tLi(:,1)), unique(length(tLp)))';
                      tLp = prdDatai.(names_fields_tL{k});
                      tp = data2plot.(names_fields_tL{k});

                      plot(tLi(:,1), tLi(:,2), '.', 'color', colors{i}, 'Markersize', 15);
                      hold on
                      plot(tp, tLp, '--', 'color', colors{i}, 'linewidth', 2);
%                       xlabel([txtData.(pets{i}).label.(names_fields_tL{k}){1}, ', ', txtData.(pets{i}).units.(names_fields_tL{k}){1}]);
%                       ylabel([txtData.(pets{i}).label.(names_fields_tL{k}){2}, ', ', txtData.(pets{i}).units.(names_fields_tL{k}){2}]); 
                      xlabel("Age (day)");
                      ylabel("Shell height (cm)"); 
                      set(gca,'Fontsize',18, 'Box', 'on')
                      print(['results_', pets{i}, '_', names_fields{j}], '-dpng')
                  end
              elseif sum(contains(names_fields, 'tL')) == 1
                  figure
                  data2plot = data.(pets{i});
                  tLi = data.(pets{i}).(names_fields{j});
                  t = linspace(min(tLi(:,1)), max(tLi(:,1)), 5000)';
                  data2plot.(names_fields{j}) = t;

    %               eval(['prdDatai = predict_', pets{i}, '(par.(pets{i}), data.(pets{i}), auxData.(pets{i}));']);
                  eval(['prdDatai = predict_', pets{i}, '(par.(pets{i}), data2plot, auxData.(pets{i}));']);
    %               tLp = unique(prdDatai.(names_fields{j}));
    %               tp = linspace(min(tLi(:,1)), max(tLi(:,1)), unique(length(tLp)))';
                  tLp = prdDatai.(names_fields{j});
                  tp = data2plot.(names_fields{j});

                  plot(tLi(:,1), tLi(:,2), '.', 'color', colors{i}, 'Markersize', 15);
                  hold on
                  plot(tp, tLp, '--', 'color', colors{i}, 'linewidth', 2);
%                   xlabel([txtData.(pets{i}).label.(names_fields{j}){1}, ', ', txtData.(pets{i}).units.(names_fields{j}){1}]);
%                   ylabel([txtData.(pets{i}).label.(names_fields{j}){2}, ', ', txtData.(pets{i}).units.(names_fields{j}){2}]);
                  xlabel("Age (day)");
                  ylabel("Shell height (cm)"); 
                  set(gca,'Fontsize',18, 'Box', 'on')
                  print(['results_', pets{i}, '_', names_fields{j}], '-dpng')
               end
          
          % weight at age
          elseif contains(names_fields{j}, 'tW') == 1
              figure
              data2plot = data.(pets{i});
%               
              tWi = data.(pets{i}).(names_fields{j});
              t = linspace(min(tWi(:,1)), max(tWi(:,1)), 5000)';
              data2plot.(names_fields{j}) = t;
              
%               eval(['prdDatai = predict_', pets{i}, '(par.(pets{i}), data.(pets{i}), auxData.(pets{i}));']);
              eval(['prdDatai = predict_', pets{i}, '(par.(pets{i}), data2plot, auxData.(pets{i}));']);
%               tWp = unique(prdDatai.(names_fields{j}));
%               tp = linspace(min(tWi(:,1)), max(tWi(:,1)), unique(length(tWp)))';
              tWp = prdDatai.(names_fields{j});
              tp = data2plot.(names_fields{j});
              
              plot(tWi(:,1), tWi(:,2), '.', 'color', colors{i}, 'Markersize', 15);
              hold on
              plot(tp, tWp, '--', 'color', colors{i}, 'linewidth', 2);
%               xlabel([txtData.(pets{i}).label.(names_fields{j}){1}, ', ', txtData.(pets{i}).units.(names_fields{j}){1}]);
%               ylabel([txtData.(pets{i}).label.(names_fields{j}){2}, ', ', txtData.(pets{i}).units.(names_fields{j}){2}]);  
              xlabel("Age (day)");
              ylabel("Tissue weight (g)"); 
              set(gca,'Fontsize',18, 'Box', 'on')
              print(['results_', pets{i}, '_', names_fields{j}], '-dpng')
              
          % weight at length   
          elseif contains(names_fields{j}, 'LW') == 1
              figure
              data2plot = data.(pets{i});
              LWi = data.(pets{i}).(names_fields{j});
              Lg = linspace(min(LWi(:,1)), max(LWi(:,1)), 5000)';
              data2plot.(names_fields{j}) = Lg;
              
%               eval(['prdDatai = predict_', pets{i}, '(par.(pets{i}), data.(pets{i}), auxData.(pets{i}));']);
              eval(['prdDatai = predict_', pets{i}, '(par.(pets{i}), data2plot, auxData.(pets{i}));']);
              LWp = prdDatai.(names_fields{j});
              Lg = data2plot.(names_fields{j});
%               LWp = unique(prdDatai.(names_fields{j}));
%               tp = linspace(min(LWi(:,1)), max(LWi(:,1)), unique(length(LWp)))';
              
              plot(LWi(:,1), LWi(:,2), '.', 'color', colors{i}, 'Markersize', 15);
              hold on
              plot(Lg, LWp, '--', 'color', colors{i}, 'linewidth', 2);
%               xlabel([txtData.(pets{i}).label.(names_fields{j}){1}, ', ', txtData.(pets{i}).units.(names_fields{j}){1}]);
%               ylabel([txtData.(pets{i}).label.(names_fields{j}){2}, ', ', txtData.(pets{i}).units.(names_fields{j}){2}]); 
              xlabel("Shell height (cm)");
              ylabel("Tissue weight (g)"); 
              set(gca,'Fontsize',18, 'Box', 'on')
              print(['results_', pets{i}, '_', names_fields{j}], '-dpng')
          end   
          
%           print -dpng results_group.png
      end 
      
  end 
  
  

  
  catch
    fprintf('Warning from custom_results_group: tis template is meant to replace the default way of presenting results\n');
    fprintf('This file requires case-specific editing, which is not yet done properly\n');
    fprintf('Use the default way of presenting results by removing this file from your current local directory\n');
    return
 end

