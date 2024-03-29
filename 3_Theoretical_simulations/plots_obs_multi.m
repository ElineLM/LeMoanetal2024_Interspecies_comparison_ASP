function plots_obs_multi(allsimu)

%------------------------------------------------------
% Objectives: plot the observable variables for several individuals simulations


% Outputs: One figure per graph, one for each observable variable
%
% calls:   allsimu, all the information and results for all individuals
% simulated


% called by:    main_multi.m

% 2024/02/06 - Eline Le Moan - based on codes from Laure Pecquerie
% Mascoet Project
%------------------------------------------------------
global nb_indiv

figure()
for i = 1:nb_indiv
    plot((allsimu(i).obs.time+allsimu(i).species.parStat.a_p)/365, allsimu(i).obs.L, 'linewidth', 2, 'color', allsimu(i).species.colorObs)
    hold on
end
xlabel('Age (year)')
ylabel('Shell height (cm)')
ylim([0 23])
xlim([0 20])
set(gca,'Fontsize',25, 'Box', 'on')
set(gcf, 'Position',  [500, 0, 800, 530])


figure()
for i = 1:nb_indiv
    plot(allsimu(i).obs.Ws, allsimu(i).obs.respi,'linewidth', 2, 'color', allsimu(i).species.colorObs)
    hold on
end
xlabel('Dry weight (g dw)')
ylabel('Respiration rate (mLO_2 h^{-1})')
ylim([0 15])
xlim([0 22])
set(gca,'Fontsize',25, 'Box', 'on')
set(gcf, 'Position',  [0, 0, 800, 530])

figure()
for i = 1:nb_indiv
    plot(allsimu(i).obs.L, allsimu(i).obs.Ws, 'linewidth', 2, 'color', allsimu(i).species.colorObs)
    hold on
end
xlabel('Shell height (cm)')
ylabel('Dry weight (g)')
ylim([0 35])
xlim([0 15])
set(gca,'Fontsize',25, 'Box', 'on')
set(gcf, 'Position',  [0, 200, 800, 530])

figure()
for i = 1:nb_indiv
    plot(allsimu(i).obs.L(end), allsimu(i).obs.meanF, 'Marker', allsimu(i).species.marker, 'MarkerSize', 10, 'linewidth', 2, 'color', allsimu(i).species.colorObs, 'MarkerFaceColor',allsimu(i).species.colorObs)
    hold on
end
xlabel('Ultimate shell height (cm)')
ylabel('Annual fecundity (# eggs)')
xlim([5 22]);
ylim([0 65000000])
ax = gca;
ax.YAxis.Exponent = 6;
set(gca,'Fontsize',25, 'Box', 'on')
set(gcf, 'Position',  [500, 200, 800, 530])
end