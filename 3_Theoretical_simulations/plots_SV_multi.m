function plots_SV_multi(allsimu)

%------------------------------------------------------
% Objectives: plot the state variables for several individuals simulations


% Outputs: One figure with 4 graphs, one for each state variable
%
% calls:   allsimu, all the information and results for all individuals
% simulated


% called by:    main_multi.m

% 2024/02/06 - Eline Le Moan - based on codes from Laure Pecquerie
% Mascoet Project
%------------------------------------------------------
global nb_indiv
    figure()

    subplot(2,2,1)
    for i = 1:nb_indiv
    plot(allsimu(i).tEVHR(:,1), allsimu(i).tEVHR(:,2), 'linewidth', 2, 'color', allsimu(i).species.colorObs)
    hold on
    end
    xlabel('Time (day)')
    ylabel('Reserve (J)')
    set(gca,'Fontsize',18, 'Box', 'on')

    subplot(2,2,2)
    for i = 1:nb_indiv
    plot(allsimu(i).tEVHR(:,1), allsimu(i).tEVHR(:,3), 'linewidth', 2, 'color', allsimu(i).species.colorObs)
    hold on
    end
    xlabel('Time (day)')
    ylabel('Structure (cm^3)')
    set(gca,'Fontsize',18, 'Box', 'on')

    subplot(2,2,3)
    for i = 1:nb_indiv
    plot(allsimu(i).tEVHR(:,1), allsimu(i).tEVHR(:,4), 'linewidth', 2, 'color', allsimu(i).species.colorObs)
    hold on
    end
    xlabel('Time (day)')
    ylabel('Maturity (J)')
    set(gca,'Fontsize',18, 'Box', 'on')

    subplot(2,2,4)
    for i = 1:nb_indiv
    plot(allsimu(i).tEVHR(:,1), allsimu(i).tEVHR(:,5), 'linewidth', 2, 'color', allsimu(i).species.colorObs)
    hold on
    end
    xlabel('Time (day)')
    ylabel('Reproduction buffer (J)')
    set(gca,'Fontsize',18, 'Box', 'on')
    set(gcf, 'Units', 'centimeters','OuterPosition', [1 1 30 22])

end