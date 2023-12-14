%% Define variables and prepare data
clear
gavr_save_folder = 'GAVR';
root_dir = 'C:/Project'; % /data/CNRL/data02/Project
components = {'P0' 'Na' 'Pa' 'Nb' 'P50' 'N1' 'P2'};
channel_cluster = {'FCz'}; % {'Frontocentral', 'Broad_Frontocentral', 'FCz'};
condition = {'11' '12' '13' '31' '32' '33' '51' '52' '53' '71' '72' '73' '91' '92' '93'};
Exp_cond = {'Quietest', 'Medium_dB', 'Loudest', 'Fastest', 'Fast',...
    'Medium_ISI', 'Slow', 'Slowest'};
dispersion_measure = 'std_err'; % std = standard deviation; std_err = standard error
color_Exp_cond = {[0 0 255], [127 127 127], [255 0 0],...
    [0 0 0], [0 0 255], [0 255 0], [255 165 0], [255 0 0]};
location_legend_P0 = 'SouthWest'; % NorthWestOutside
location_legend_Na = 'SouthWest'; % NorthWestOutside
location_legend_Pa = 'SouthWest'; % NorthWestOutside
location_legend_Nb = 'SouthWest'; % NorthWestOutside
location_legend_P50 = 'NorthWest'; % NorthWestOutside
location_legend_N1 = 'SouthWest'; % NorthWestOutside
location_legend_P2 = 'NorthWest';
legend_fontsize = 9; % normally is 6
marker_size = 6;
axes_y_P2 = [1 7]; axes_y_N1 = [-3.5 0]; axes_y_P50 = [0.1 1.3]; axes_y_Nb = [-0.4 0.2];
axes_y_Pa = [0.1 0.7]; axes_y_Na = [-0.6 -0.1]; axes_y_P0 = [-0.1 0.5];
single_or_combined = 2; % 1 = separated plots; 2 = combined into a 2x2 figure

% Prepare data for all components
for chan = 1:length(channel_cluster)
    for cm = 1:length(components)
        load ([root_dir '/Statistics/Graph_components/' gavr_save_folder '/' components{cm} '_' channel_cluster{chan} '.mat']);
        for c = 1:length(condition)
            eval(['average_' components{cm} '_' condition{c} ' = mean(' components{cm} '_' condition{c} ',2);'])
            if strcmp(dispersion_measure,'std')
                eval(['disp_' components{cm} '_' condition{c} ' = std(' components{cm} '_' condition{c} ',0,2);'])
            elseif strcmp(dispersion_measure,'std_err') % standard error
                eval(['disp_' components{cm} '_' condition{c} ' = (std(' components{cm} '_' condition{c} ',0,2))/sqrt(size(' components{cm} '_' condition{c} ',2));'])
            end
        end
    end
end

% Prepare data for averaged conditions
for chan = 1:length(channel_cluster)
    for cm = 1:length(components)
        load ([root_dir '/Statistics/Graph_components/' gavr_save_folder '/' components{cm} '_' channel_cluster{chan} '.mat']);
        for ex = 1:length(Exp_cond)
            eval(['average_' components{cm} '_' Exp_cond{ex} ' = mean(' components{cm} '_' Exp_cond{ex} ',2);'])
            if strcmp(dispersion_measure,'std')
                eval(['disp_' components{cm} '_' Exp_cond{ex} ' = std(' components{cm} '_' Exp_cond{ex} ',0,2);'])
            elseif strcmp(dispersion_measure,'std_err') % standard error
                eval(['disp_' components{cm} '_' Exp_cond{ex} ' = (std(' components{cm} '_' Exp_cond{ex} ',0,2))/sqrt(size(' components{cm} '_' Exp_cond{ex} ',2));'])
            end
        end
    end
end

%% P2
if single_or_combined == 1
    figure;
else
    figure('units','normalized','outerposition',[0 0 1 1]);
end
% Divided by Intensity
Grafic_P2 = ...
     [average_P2_11	average_P2_31 average_P2_51 average_P2_71 average_P2_91;...
     disp_P2_11	disp_P2_31 disp_P2_51 disp_P2_71 disp_P2_91; ...
     average_P2_12	average_P2_32 average_P2_52 average_P2_72 average_P2_92; ... 
     disp_P2_12	disp_P2_32 disp_P2_52 disp_P2_72 disp_P2_92;...
     average_P2_13	average_P2_33 average_P2_53 average_P2_73 average_P2_93; ... 
     disp_P2_13	disp_P2_33 disp_P2_53 disp_P2_73 disp_P2_93];
 
 hold on; 
 if single_or_combined == 2
    subplot(2,2,1);
 end
 plot(Grafic_P2(1,:),'color', color_Exp_cond{1}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P2(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{1}/256,'MarkerFaceColor',color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P2(1,:),Grafic_P2(2,:),'color', color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P2(3,:),'color', color_Exp_cond{2}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P2(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{2}/256,'MarkerFaceColor',color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P2(3,:),Grafic_P2(4,:),'color', color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P2(5,:),'color', color_Exp_cond{3}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P2(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{3}/256,'MarkerFaceColor',color_Exp_cond{3}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P2(5,:),Grafic_P2(6,:),'color', color_Exp_cond{3}/256,'HandleVisibility','off')
 lgd = legend('65dB','75dB','85dB','location', location_legend_P2);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes P2')
set(gca,'XLim',[0.5 5.5]); 
set(gca,'YLim',axes_y_P2); 
%xticklabels({[],'0.25-0.5',[],'0.5-1',[],'1-2',[],'2-4',[],'4-8',[]})
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3 4 5]; 
xticklabels({'0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Divided by ISI
Grafic_P2 = ...
     [average_P2_11	average_P2_12 average_P2_13;...
     disp_P2_11	disp_P2_12 disp_P2_13; ...
     average_P2_31	average_P2_32 average_P2_33;...
     disp_P2_31	disp_P2_32 disp_P2_33; ...
     average_P2_51	average_P2_52 average_P2_53;...
     disp_P2_51	disp_P2_52 disp_P2_53; ...
     average_P2_71	average_P2_72 average_P2_73;...
     disp_P2_71	disp_P2_72 disp_P2_73; ...
     average_P2_91	average_P2_92 average_P2_93;...
     disp_P2_91	disp_P2_92 disp_P2_93];
 
 if single_or_combined == 1
    figure;
 end
 hold on;
 if single_or_combined == 2
    subplot(2,2,2);
 end
 plot(Grafic_P2(1,:),'color', color_Exp_cond{4}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P2(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{4}/256,'MarkerFaceColor',color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P2(1,:),Grafic_P2(2,:),'color', color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P2(3,:),'color', color_Exp_cond{5}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P2(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{5}/256,'MarkerFaceColor',color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P2(3,:),Grafic_P2(4,:),'color', color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P2(5,:),'color', color_Exp_cond{6}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P2(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{6}/256,'MarkerFaceColor',color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P2(5,:),Grafic_P2(6,:),'color', color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P2(7,:),'color', color_Exp_cond{7}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P2(7,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{7}/256,'MarkerFaceColor',color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P2(7,:),Grafic_P2(8,:),'color', color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P2(9,:),'color', color_Exp_cond{8}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P2(9,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{8}/256,'MarkerFaceColor',color_Exp_cond{8}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P2(9,:),Grafic_P2(10,:),'color', color_Exp_cond{8}/256,'HandleVisibility','off')
 
lgd = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','location', location_legend_P2);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes P2')
set(gca,'XLim',[0.75 3.25]); 
set(gca,'YLim',axes_y_P2); 
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3]; 
xticklabels({'65dB','75dB','85dB'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Averaged across ISI
% Data to be plotted as a bar graph
model_series = [average_P2_Quietest; average_P2_Medium_dB; average_P2_Loudest];
%Data to be plotted as the error bars
model_error = [disp_P2_Quietest; disp_P2_Medium_dB; disp_P2_Loudest]; 

if single_or_combined == 1
    figure;
else
    subplot(2,2,3);
end

% Creating axes and the bar graph
hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('65dB','75dB','85dB','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
title('Amplitudes P2 averaged across ISI');



% Averaged across dB
model_series = [average_P2_Fastest; average_P2_Fast; average_P2_Medium_ISI; average_P2_Slow; average_P2_Slowest];
%Data to be plotted as the error bars
model_error = [disp_P2_Fastest; disp_P2_Fast; disp_P2_Medium_ISI; disp_P2_Slow; disp_P2_Slowest]; 

% Creating axes and the bar graph
if single_or_combined == 1
    figure;
else
    subplot(2,2,4);
end

hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i+3}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
errorbar(x(i), model_series(1,1), model_error(1,1), 'r', 'linestyle', 'none','HandleVisibility','off');
title('Amplitudes P2 averaged across dB');

%% N1

if single_or_combined == 1
    figure;
else
    figure('units','normalized','outerposition',[0 0 1 1]);
end
% Divided by Intensity
Grafic_N1 = ...
     [average_N1_11	average_N1_31 average_N1_51 average_N1_71 average_N1_91;...
     disp_N1_11	disp_N1_31 disp_N1_51 disp_N1_71 disp_N1_91; ...
     average_N1_12	average_N1_32 average_N1_52 average_N1_72 average_N1_92; ... 
     disp_N1_12	disp_N1_32 disp_N1_52 disp_N1_72 disp_N1_92;...
     average_N1_13	average_N1_33 average_N1_53 average_N1_73 average_N1_93; ... 
     disp_N1_13	disp_N1_33 disp_N1_53 disp_N1_73 disp_N1_93];
 
 hold on; 
 if single_or_combined == 2
    subplot(2,2,1);
 end
 plot(Grafic_N1(1,:),'color', color_Exp_cond{1}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_N1(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{1}/256,'MarkerFaceColor',color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_N1(1,:),Grafic_N1(2,:),'color', color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; plot(Grafic_N1(3,:),'color', color_Exp_cond{2}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_N1(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{2}/256,'MarkerFaceColor',color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_N1(3,:),Grafic_N1(4,:),'color', color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; plot(Grafic_N1(5,:),'color', color_Exp_cond{3}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_N1(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{3}/256,'MarkerFaceColor',color_Exp_cond{3}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_N1(5,:),Grafic_N1(6,:),'color', color_Exp_cond{3}/256,'HandleVisibility','off')
 lgd = legend('65dB','75dB','85dB','location', location_legend_N1);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes N1')
set(gca,'XLim',[0.5 5.5]); 
set(gca,'YLim',axes_y_N1); 
%xticklabels({[],'0.25-0.5',[],'0.5-1',[],'1-2',[],'2-4',[],'4-8',[]})
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3 4 5]; 
xticklabels({'0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Divided by ISI
Grafic_N1 = ...
     [average_N1_11	average_N1_12 average_N1_13;...
     disp_N1_11	disp_N1_12 disp_N1_13; ...
     average_N1_31	average_N1_32 average_N1_33;...
     disp_N1_31	disp_N1_32 disp_N1_33; ...
     average_N1_51	average_N1_52 average_N1_53;...
     disp_N1_51	disp_N1_52 disp_N1_53; ...
     average_N1_71	average_N1_72 average_N1_73;...
     disp_N1_71	disp_N1_72 disp_N1_73; ...
     average_N1_91	average_N1_92 average_N1_93;...
     disp_N1_91	disp_N1_92 disp_N1_93];
 
 if single_or_combined == 1
    figure;
 end
 hold on;
 if single_or_combined == 2
    subplot(2,2,2);
 end
 plot(Grafic_N1(1,:),'color', color_Exp_cond{4}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_N1(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{4}/256,'MarkerFaceColor',color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_N1(1,:),Grafic_N1(2,:),'color', color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; plot(Grafic_N1(3,:),'color', color_Exp_cond{5}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_N1(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{5}/256,'MarkerFaceColor',color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_N1(3,:),Grafic_N1(4,:),'color', color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; plot(Grafic_N1(5,:),'color', color_Exp_cond{6}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_N1(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{6}/256,'MarkerFaceColor',color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_N1(5,:),Grafic_N1(6,:),'color', color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; plot(Grafic_N1(7,:),'color', color_Exp_cond{7}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_N1(7,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{7}/256,'MarkerFaceColor',color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_N1(7,:),Grafic_N1(8,:),'color', color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; plot(Grafic_N1(9,:),'color', color_Exp_cond{8}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_N1(9,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{8}/256,'MarkerFaceColor',color_Exp_cond{8}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_N1(9,:),Grafic_N1(10,:),'color', color_Exp_cond{8}/256,'HandleVisibility','off')
 
lgd = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','location', location_legend_N1);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes N1')
set(gca,'XLim',[0.75 3.25]); 
set(gca,'YLim',axes_y_N1); 
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3]; 
xticklabels({'65dB','75dB','85dB'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Averaged across ISI
% Data to be plotted as a bar graph
model_series = [average_N1_Quietest; average_N1_Medium_dB; average_N1_Loudest];
%Data to be plotted as the error bars
model_error = [disp_N1_Quietest; disp_N1_Medium_dB; disp_N1_Loudest]; 

if single_or_combined == 1
    figure;
else
    subplot(2,2,3);
end

% Creating axes and the bar graph
hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('65dB','75dB','85dB','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
title('Amplitudes N1 averaged across ISI');



% Averaged across dB
model_series = [average_N1_Fastest; average_N1_Fast; average_N1_Medium_ISI; average_N1_Slow; average_N1_Slowest];
%Data to be plotted as the error bars
model_error = [disp_N1_Fastest; disp_N1_Fast; disp_N1_Medium_ISI; disp_N1_Slow; disp_N1_Slowest]; 

% Creating axes and the bar graph
if single_or_combined == 1
    figure;
else
    subplot(2,2,4);
end

hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i+3}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
errorbar(x(i), model_series(1,1), model_error(1,1), 'r', 'linestyle', 'none','HandleVisibility','off');
title('Amplitudes N1 averaged across dB');

%% P50
if single_or_combined == 1
    figure;
else
    figure('units','normalized','outerposition',[0 0 1 1]);
end
% Divided by Intensity
Grafic_P50 = ...
     [average_P50_11	average_P50_31 average_P50_51 average_P50_71 average_P50_91;...
     disp_P50_11	disp_P50_31 disp_P50_51 disp_P50_71 disp_P50_91; ...
     average_P50_12	average_P50_32 average_P50_52 average_P50_72 average_P50_92; ... 
     disp_P50_12	disp_P50_32 disp_P50_52 disp_P50_72 disp_P50_92;...
     average_P50_13	average_P50_33 average_P50_53 average_P50_73 average_P50_93; ... 
     disp_P50_13	disp_P50_33 disp_P50_53 disp_P50_73 disp_P50_93];
 
 hold on; 
 if single_or_combined == 2
    subplot(2,2,1);
 end
 plot(Grafic_P50(1,:),'color', color_Exp_cond{1}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P50(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{1}/256,'MarkerFaceColor',color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P50(1,:),Grafic_P50(2,:),'color', color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P50(3,:),'color', color_Exp_cond{2}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P50(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{2}/256,'MarkerFaceColor',color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P50(3,:),Grafic_P50(4,:),'color', color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P50(5,:),'color', color_Exp_cond{3}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P50(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{3}/256,'MarkerFaceColor',color_Exp_cond{3}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P50(5,:),Grafic_P50(6,:),'color', color_Exp_cond{3}/256,'HandleVisibility','off')
 lgd = legend('65dB','75dB','85dB','location', location_legend_P50);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes P50')
set(gca,'XLim',[0.5 5.5]); 
set(gca,'YLim',axes_y_P50); 
%xticklabels({[],'0.25-0.5',[],'0.5-1',[],'1-2',[],'2-4',[],'4-8',[]})
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3 4 5]; 
xticklabels({'0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Divided by ISI
Grafic_P50 = ...
     [average_P50_11	average_P50_12 average_P50_13;...
     disp_P50_11	disp_P50_12 disp_P50_13; ...
     average_P50_31	average_P50_32 average_P50_33;...
     disp_P50_31	disp_P50_32 disp_P50_33; ...
     average_P50_51	average_P50_52 average_P50_53;...
     disp_P50_51	disp_P50_52 disp_P50_53; ...
     average_P50_71	average_P50_72 average_P50_73;...
     disp_P50_71	disp_P50_72 disp_P50_73; ...
     average_P50_91	average_P50_92 average_P50_93;...
     disp_P50_91	disp_P50_92 disp_P50_93];
 
 if single_or_combined == 1
    figure;
 end
 hold on;
 if single_or_combined == 2
    subplot(2,2,2);
 end
 plot(Grafic_P50(1,:),'color', color_Exp_cond{4}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P50(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{4}/256,'MarkerFaceColor',color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P50(1,:),Grafic_P50(2,:),'color', color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P50(3,:),'color', color_Exp_cond{5}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P50(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{5}/256,'MarkerFaceColor',color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P50(3,:),Grafic_P50(4,:),'color', color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P50(5,:),'color', color_Exp_cond{6}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P50(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{6}/256,'MarkerFaceColor',color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P50(5,:),Grafic_P50(6,:),'color', color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P50(7,:),'color', color_Exp_cond{7}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P50(7,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{7}/256,'MarkerFaceColor',color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P50(7,:),Grafic_P50(8,:),'color', color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P50(9,:),'color', color_Exp_cond{8}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P50(9,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{8}/256,'MarkerFaceColor',color_Exp_cond{8}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P50(9,:),Grafic_P50(10,:),'color', color_Exp_cond{8}/256,'HandleVisibility','off')
 
lgd = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','location', location_legend_P50);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes P50')
set(gca,'XLim',[0.75 3.25]); 
set(gca,'YLim',axes_y_P50); 
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3]; 
xticklabels({'65dB','75dB','85dB'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Averaged across ISI
% Data to be plotted as a bar graph
model_series = [average_P50_Quietest; average_P50_Medium_dB; average_P50_Loudest];
%Data to be plotted as the error bars
model_error = [disp_P50_Quietest; disp_P50_Medium_dB; disp_P50_Loudest]; 

if single_or_combined == 1
    figure;
else
    subplot(2,2,3);
end

% Creating axes and the bar graph
hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('65dB','75dB','85dB','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
title('Amplitudes P50 averaged across ISI');



% Averaged across dB
model_series = [average_P50_Fastest; average_P50_Fast; average_P50_Medium_ISI; average_P50_Slow; average_P50_Slowest];
%Data to be plotted as the error bars
model_error = [disp_P50_Fastest; disp_P50_Fast; disp_P50_Medium_ISI; disp_P50_Slow; disp_P50_Slowest]; 

% Creating axes and the bar graph
if single_or_combined == 1
    figure;
else
    subplot(2,2,4);
end

hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i+3}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
errorbar(x(i), model_series(1,1), model_error(1,1), 'r', 'linestyle', 'none','HandleVisibility','off');
title('Amplitudes P50 averaged across dB');

%% Nb
 if single_or_combined == 1
    figure;
else
    figure('units','normalized','outerposition',[0 0 1 1]);
end
% Divided by Intensity
Grafic_Nb = ...
     [average_Nb_11	average_Nb_31 average_Nb_51 average_Nb_71 average_Nb_91;...
     disp_Nb_11	disp_Nb_31 disp_Nb_51 disp_Nb_71 disp_Nb_91; ...
     average_Nb_12	average_Nb_32 average_Nb_52 average_Nb_72 average_Nb_92; ... 
     disp_Nb_12	disp_Nb_32 disp_Nb_52 disp_Nb_72 disp_Nb_92;...
     average_Nb_13	average_Nb_33 average_Nb_53 average_Nb_73 average_Nb_93; ... 
     disp_Nb_13	disp_Nb_33 disp_Nb_53 disp_Nb_73 disp_Nb_93];
 
 hold on; 
 if single_or_combined == 2
    subplot(2,2,1);
 end
 plot(Grafic_Nb(1,:),'color', color_Exp_cond{1}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Nb(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{1}/256,'MarkerFaceColor',color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Nb(1,:),Grafic_Nb(2,:),'color', color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Nb(3,:),'color', color_Exp_cond{2}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Nb(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{2}/256,'MarkerFaceColor',color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Nb(3,:),Grafic_Nb(4,:),'color', color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Nb(5,:),'color', color_Exp_cond{3}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Nb(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{3}/256,'MarkerFaceColor',color_Exp_cond{3}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Nb(5,:),Grafic_Nb(6,:),'color', color_Exp_cond{3}/256,'HandleVisibility','off')
 lgd = legend('65dB','75dB','85dB','location', location_legend_Nb);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes Nb')
set(gca,'XLim',[0.5 5.5]); 
set(gca,'YLim',axes_y_Nb); 
%xticklabels({[],'0.25-0.5',[],'0.5-1',[],'1-2',[],'2-4',[],'4-8',[]})
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3 4 5]; 
xticklabels({'0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Divided by ISI
Grafic_Nb = ...
     [average_Nb_11	average_Nb_12 average_Nb_13;...
     disp_Nb_11	disp_Nb_12 disp_Nb_13; ...
     average_Nb_31	average_Nb_32 average_Nb_33;...
     disp_Nb_31	disp_Nb_32 disp_Nb_33; ...
     average_Nb_51	average_Nb_52 average_Nb_53;...
     disp_Nb_51	disp_Nb_52 disp_Nb_53; ...
     average_Nb_71	average_Nb_72 average_Nb_73;...
     disp_Nb_71	disp_Nb_72 disp_Nb_73; ...
     average_Nb_91	average_Nb_92 average_Nb_93;...
     disp_Nb_91	disp_Nb_92 disp_Nb_93];
 
 if single_or_combined == 1
    figure;
 end
 hold on;
 if single_or_combined == 2
    subplot(2,2,2);
 end
 plot(Grafic_Nb(1,:),'color', color_Exp_cond{4}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Nb(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{4}/256,'MarkerFaceColor',color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Nb(1,:),Grafic_Nb(2,:),'color', color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Nb(3,:),'color', color_Exp_cond{5}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Nb(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{5}/256,'MarkerFaceColor',color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Nb(3,:),Grafic_Nb(4,:),'color', color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Nb(5,:),'color', color_Exp_cond{6}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Nb(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{6}/256,'MarkerFaceColor',color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Nb(5,:),Grafic_Nb(6,:),'color', color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Nb(7,:),'color', color_Exp_cond{7}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Nb(7,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{7}/256,'MarkerFaceColor',color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Nb(7,:),Grafic_Nb(8,:),'color', color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Nb(9,:),'color', color_Exp_cond{8}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Nb(9,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{8}/256,'MarkerFaceColor',color_Exp_cond{8}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Nb(9,:),Grafic_Nb(10,:),'color', color_Exp_cond{8}/256,'HandleVisibility','off')
 
lgd = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','location', location_legend_Nb);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes Nb')
set(gca,'XLim',[0.75 3.25]); 
set(gca,'YLim',axes_y_Nb); 
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3]; 
xticklabels({'65dB','75dB','85dB'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Averaged across ISI
% Data to be plotted as a bar graph
model_series = [average_Nb_Quietest; average_Nb_Medium_dB; average_Nb_Loudest];
%Data to be plotted as the error bars
model_error = [disp_Nb_Quietest; disp_Nb_Medium_dB; disp_Nb_Loudest]; 

if single_or_combined == 1
    figure;
else
    subplot(2,2,3);
end

% Creating axes and the bar graph
hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('65dB','75dB','85dB','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
title('Amplitudes Nb averaged across ISI');



% Averaged across dB
model_series = [average_Nb_Fastest; average_Nb_Fast; average_Nb_Medium_ISI; average_Nb_Slow; average_Nb_Slowest];
%Data to be plotted as the error bars
model_error = [disp_Nb_Fastest; disp_Nb_Fast; disp_Nb_Medium_ISI; disp_Nb_Slow; disp_Nb_Slowest]; 

% Creating axes and the bar graph
if single_or_combined == 1
    figure;
else
    subplot(2,2,4);
end

hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i+3}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
errorbar(x(i), model_series(1,1), model_error(1,1), 'r', 'linestyle', 'none','HandleVisibility','off');
title('Amplitudes Nb averaged across dB');

%% Pa
 if single_or_combined == 1
    figure;
else
    figure('units','normalized','outerposition',[0 0 1 1]);
end
% Divided by Intensity
Grafic_Pa = ...
     [average_Pa_11	average_Pa_31 average_Pa_51 average_Pa_71 average_Pa_91;...
     disp_Pa_11	disp_Pa_31 disp_Pa_51 disp_Pa_71 disp_Pa_91; ...
     average_Pa_12	average_Pa_32 average_Pa_52 average_Pa_72 average_Pa_92; ... 
     disp_Pa_12	disp_Pa_32 disp_Pa_52 disp_Pa_72 disp_Pa_92;...
     average_Pa_13	average_Pa_33 average_Pa_53 average_Pa_73 average_Pa_93; ... 
     disp_Pa_13	disp_Pa_33 disp_Pa_53 disp_Pa_73 disp_Pa_93];
 
 hold on; 
 if single_or_combined == 2
    subplot(2,2,1);
 end
 plot(Grafic_Pa(1,:),'color', color_Exp_cond{1}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Pa(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{1}/256,'MarkerFaceColor',color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Pa(1,:),Grafic_Pa(2,:),'color', color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Pa(3,:),'color', color_Exp_cond{2}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Pa(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{2}/256,'MarkerFaceColor',color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Pa(3,:),Grafic_Pa(4,:),'color', color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Pa(5,:),'color', color_Exp_cond{3}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Pa(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{3}/256,'MarkerFaceColor',color_Exp_cond{3}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Pa(5,:),Grafic_Pa(6,:),'color', color_Exp_cond{3}/256,'HandleVisibility','off')
 lgd = legend('65dB','75dB','85dB','location', location_legend_Pa);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes Pa')
set(gca,'XLim',[0.5 5.5]); 
set(gca,'YLim',axes_y_Pa); 
%xticklabels({[],'0.25-0.5',[],'0.5-1',[],'1-2',[],'2-4',[],'4-8',[]})
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3 4 5]; 
xticklabels({'0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Divided by ISI
Grafic_Pa = ...
     [average_Pa_11	average_Pa_12 average_Pa_13;...
     disp_Pa_11	disp_Pa_12 disp_Pa_13; ...
     average_Pa_31	average_Pa_32 average_Pa_33;...
     disp_Pa_31	disp_Pa_32 disp_Pa_33; ...
     average_Pa_51	average_Pa_52 average_Pa_53;...
     disp_Pa_51	disp_Pa_52 disp_Pa_53; ...
     average_Pa_71	average_Pa_72 average_Pa_73;...
     disp_Pa_71	disp_Pa_72 disp_Pa_73; ...
     average_Pa_91	average_Pa_92 average_Pa_93;...
     disp_Pa_91	disp_Pa_92 disp_Pa_93];
 
 if single_or_combined == 1
    figure;
 end
 hold on;
 if single_or_combined == 2
    subplot(2,2,2);
 end
 plot(Grafic_Pa(1,:),'color', color_Exp_cond{4}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Pa(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{4}/256,'MarkerFaceColor',color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Pa(1,:),Grafic_Pa(2,:),'color', color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Pa(3,:),'color', color_Exp_cond{5}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Pa(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{5}/256,'MarkerFaceColor',color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Pa(3,:),Grafic_Pa(4,:),'color', color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Pa(5,:),'color', color_Exp_cond{6}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Pa(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{6}/256,'MarkerFaceColor',color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Pa(5,:),Grafic_Pa(6,:),'color', color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Pa(7,:),'color', color_Exp_cond{7}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Pa(7,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{7}/256,'MarkerFaceColor',color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Pa(7,:),Grafic_Pa(8,:),'color', color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Pa(9,:),'color', color_Exp_cond{8}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Pa(9,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{8}/256,'MarkerFaceColor',color_Exp_cond{8}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Pa(9,:),Grafic_Pa(10,:),'color', color_Exp_cond{8}/256,'HandleVisibility','off')
 
lgd = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','location', location_legend_Pa);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes Pa')
set(gca,'XLim',[0.75 3.25]); 
set(gca,'YLim',axes_y_Pa); 
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3]; 
xticklabels({'65dB','75dB','85dB'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Averaged across ISI
% Data to be plotted as a bar graph
model_series = [average_Pa_Quietest; average_Pa_Medium_dB; average_Pa_Loudest];
%Data to be plotted as the error bars
model_error = [disp_Pa_Quietest; disp_Pa_Medium_dB; disp_Pa_Loudest]; 

if single_or_combined == 1
    figure;
else
    subplot(2,2,3);
end

% Creating axes and the bar graph
hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('65dB','75dB','85dB','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
title('Amplitudes Pa averaged across ISI');



% Averaged across dB
model_series = [average_Pa_Fastest; average_Pa_Fast; average_Pa_Medium_ISI; average_Pa_Slow; average_Pa_Slowest];
%Data to be plotted as the error bars
model_error = [disp_Pa_Fastest; disp_Pa_Fast; disp_Pa_Medium_ISI; disp_Pa_Slow; disp_Pa_Slowest]; 

% Creating axes and the bar graph
if single_or_combined == 1
    figure;
else
    subplot(2,2,4);
end

hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i+3}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
errorbar(x(i), model_series(1,1), model_error(1,1), 'r', 'linestyle', 'none','HandleVisibility','off');
title('Amplitudes Pa averaged across dB');

 %% Na
 if single_or_combined == 1
    figure;
else
    figure('units','normalized','outerposition',[0 0 1 1]);
end
% Divided by Intensity
Grafic_Na = ...
     [average_Na_11	average_Na_31 average_Na_51 average_Na_71 average_Na_91;...
     disp_Na_11	disp_Na_31 disp_Na_51 disp_Na_71 disp_Na_91; ...
     average_Na_12	average_Na_32 average_Na_52 average_Na_72 average_Na_92; ... 
     disp_Na_12	disp_Na_32 disp_Na_52 disp_Na_72 disp_Na_92;...
     average_Na_13	average_Na_33 average_Na_53 average_Na_73 average_Na_93; ... 
     disp_Na_13	disp_Na_33 disp_Na_53 disp_Na_73 disp_Na_93];
 
 hold on; 
 if single_or_combined == 2
    subplot(2,2,1);
 end
 plot(Grafic_Na(1,:),'color', color_Exp_cond{1}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Na(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{1}/256,'MarkerFaceColor',color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Na(1,:),Grafic_Na(2,:),'color', color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Na(3,:),'color', color_Exp_cond{2}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Na(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{2}/256,'MarkerFaceColor',color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Na(3,:),Grafic_Na(4,:),'color', color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Na(5,:),'color', color_Exp_cond{3}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Na(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{3}/256,'MarkerFaceColor',color_Exp_cond{3}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_Na(5,:),Grafic_Na(6,:),'color', color_Exp_cond{3}/256,'HandleVisibility','off')
 lgd = legend('65dB','75dB','85dB','location', location_legend_Na);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes Na')
set(gca,'XLim',[0.5 5.5]); 
set(gca,'YLim',axes_y_Na); 
%xticklabels({[],'0.25-0.5',[],'0.5-1',[],'1-2',[],'2-4',[],'4-8',[]})
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3 4 5]; 
xticklabels({'0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Divided by ISI
Grafic_Na = ...
     [average_Na_11	average_Na_12 average_Na_13;...
     disp_Na_11	disp_Na_12 disp_Na_13; ...
     average_Na_31	average_Na_32 average_Na_33;...
     disp_Na_31	disp_Na_32 disp_Na_33; ...
     average_Na_51	average_Na_52 average_Na_53;...
     disp_Na_51	disp_Na_52 disp_Na_53; ...
     average_Na_71	average_Na_72 average_Na_73;...
     disp_Na_71	disp_Na_72 disp_Na_73; ...
     average_Na_91	average_Na_92 average_Na_93;...
     disp_Na_91	disp_Na_92 disp_Na_93];
 
 if single_or_combined == 1
    figure;
 end
 hold on;
 if single_or_combined == 2
    subplot(2,2,2);
 end
 plot(Grafic_Na(1,:),'color', color_Exp_cond{4}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Na(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{4}/256,'MarkerFaceColor',color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Na(1,:),Grafic_Na(2,:),'color', color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Na(3,:),'color', color_Exp_cond{5}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Na(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{5}/256,'MarkerFaceColor',color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Na(3,:),Grafic_Na(4,:),'color', color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Na(5,:),'color', color_Exp_cond{6}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Na(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{6}/256,'MarkerFaceColor',color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Na(5,:),Grafic_Na(6,:),'color', color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Na(7,:),'color', color_Exp_cond{7}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Na(7,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{7}/256,'MarkerFaceColor',color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Na(7,:),Grafic_Na(8,:),'color', color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; plot(Grafic_Na(9,:),'color', color_Exp_cond{8}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_Na(9,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{8}/256,'MarkerFaceColor',color_Exp_cond{8}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_Na(9,:),Grafic_Na(10,:),'color', color_Exp_cond{8}/256,'HandleVisibility','off')
 
lgd = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','location', location_legend_Na);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes Na')
set(gca,'XLim',[0.75 3.25]); 
set(gca,'YLim',axes_y_Na); 
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3]; 
xticklabels({'65dB','75dB','85dB'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Averaged across ISI
% Data to be plotted as a bar graph
model_series = [average_Na_Quietest; average_Na_Medium_dB; average_Na_Loudest];
%Data to be plotted as the error bars
model_error = [disp_Na_Quietest; disp_Na_Medium_dB; disp_Na_Loudest]; 

if single_or_combined == 1
    figure;
else
    subplot(2,2,3);
end

% Creating axes and the bar graph
hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('65dB','75dB','85dB','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
title('Amplitudes Na averaged across ISI');



% Averaged across dB
model_series = [average_Na_Fastest; average_Na_Fast; average_Na_Medium_ISI; average_Na_Slow; average_Na_Slowest];
%Data to be plotted as the error bars
model_error = [disp_Na_Fastest; disp_Na_Fast; disp_Na_Medium_ISI; disp_Na_Slow; disp_Na_Slowest]; 

% Creating axes and the bar graph
if single_or_combined == 1
    figure;
else
    subplot(2,2,4);
end

hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i+3}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
errorbar(x(i), model_series(1,1), model_error(1,1), 'r', 'linestyle', 'none','HandleVisibility','off');
title('Amplitudes Na averaged across dB');

 %% P0
if single_or_combined == 1
    figure;
else
    figure('units','normalized','outerposition',[0 0 1 1]);
end
% Divided by Intensity
Grafic_P0 = ...
     [average_P0_11	average_P0_31 average_P0_51 average_P0_71 average_P0_91;...
     disp_P0_11	disp_P0_31 disp_P0_51 disp_P0_71 disp_P0_91; ...
     average_P0_12	average_P0_32 average_P0_52 average_P0_72 average_P0_92; ... 
     disp_P0_12	disp_P0_32 disp_P0_52 disp_P0_72 disp_P0_92;...
     average_P0_13	average_P0_33 average_P0_53 average_P0_73 average_P0_93; ... 
     disp_P0_13	disp_P0_33 disp_P0_53 disp_P0_73 disp_P0_93];
 
 hold on; 
 if single_or_combined == 2
    subplot(2,2,1);
 end
 plot(Grafic_P0(1,:),'color', color_Exp_cond{1}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P0(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{1}/256,'MarkerFaceColor',color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P0(1,:),Grafic_P0(2,:),'color', color_Exp_cond{1}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P0(3,:),'color', color_Exp_cond{2}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P0(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{2}/256,'MarkerFaceColor',color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P0(3,:),Grafic_P0(4,:),'color', color_Exp_cond{2}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P0(5,:),'color', color_Exp_cond{3}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P0(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{3}/256,'MarkerFaceColor',color_Exp_cond{3}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3 4 5],Grafic_P0(5,:),Grafic_P0(6,:),'color', color_Exp_cond{3}/256,'HandleVisibility','off')
 lgd = legend('65dB','75dB','85dB','location', location_legend_P0);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes P0')
set(gca,'XLim',[0.5 5.5]); 
set(gca,'YLim',axes_y_P0); 
%xticklabels({[],'0.25-0.5',[],'0.5-1',[],'1-2',[],'2-4',[],'4-8',[]})
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3 4 5]; 
xticklabels({'0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Divided by ISI
Grafic_P0 = ...
     [average_P0_11	average_P0_12 average_P0_13;...
     disp_P0_11	disp_P0_12 disp_P0_13; ...
     average_P0_31	average_P0_32 average_P0_33;...
     disp_P0_31	disp_P0_32 disp_P0_33; ...
     average_P0_51	average_P0_52 average_P0_53;...
     disp_P0_51	disp_P0_52 disp_P0_53; ...
     average_P0_71	average_P0_72 average_P0_73;...
     disp_P0_71	disp_P0_72 disp_P0_73; ...
     average_P0_91	average_P0_92 average_P0_93;...
     disp_P0_91	disp_P0_92 disp_P0_93];
 
 if single_or_combined == 1
    figure;
 end
 hold on;
 if single_or_combined == 2
    subplot(2,2,2);
 end
 plot(Grafic_P0(1,:),'color', color_Exp_cond{4}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P0(1,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{4}/256,'MarkerFaceColor',color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P0(1,:),Grafic_P0(2,:),'color', color_Exp_cond{4}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P0(3,:),'color', color_Exp_cond{5}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P0(3,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{5}/256,'MarkerFaceColor',color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P0(3,:),Grafic_P0(4,:),'color', color_Exp_cond{5}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P0(5,:),'color', color_Exp_cond{6}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P0(5,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{6}/256,'MarkerFaceColor',color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P0(5,:),Grafic_P0(6,:),'color', color_Exp_cond{6}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P0(7,:),'color', color_Exp_cond{7}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P0(7,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{7}/256,'MarkerFaceColor',color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P0(7,:),Grafic_P0(8,:),'color', color_Exp_cond{7}/256,'HandleVisibility','off')
 hold on; plot(Grafic_P0(9,:),'color', color_Exp_cond{8}/256, 'LineWidth', 1.5);
 hold on; plot(Grafic_P0(9,:),'-s','MarkerSize',marker_size,'MarkerEdgeColor',color_Exp_cond{8}/256,'MarkerFaceColor',color_Exp_cond{8}/256,'HandleVisibility','off')
 hold on; errorbar([1 2 3],Grafic_P0(9,:),Grafic_P0(10,:),'color', color_Exp_cond{8}/256,'HandleVisibility','off')
 
lgd = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','location', location_legend_P0);
lgd.FontSize = legend_fontsize;
lgd.EdgeColor = [1 1 1];
lgd.ItemTokenSize = [10 18];
title('Amplitudes P0')
set(gca,'XLim',[0.75 3.25]); 
set(gca,'YLim',axes_y_P0); 
hAx = ancestor(gca, 'axes');
xt=hAx.XTick;
hAx.XTick=[1 2 3]; 
xticklabels({'65dB','75dB','85dB'})
ylabel('Amplitude (µV)'); %µV2/Hz

% Averaged across ISI
% Data to be plotted as a bar graph
model_series = [average_P0_Quietest; average_P0_Medium_dB; average_P0_Loudest];
%Data to be plotted as the error bars
model_error = [disp_P0_Quietest; disp_P0_Medium_dB; disp_P0_Loudest]; 

if single_or_combined == 1
    figure;
else
    subplot(2,2,3);
end

% Creating axes and the bar graph
hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('65dB','75dB','85dB','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
title('Amplitudes P0 averaged across ISI');



% Averaged across dB
model_series = [average_P0_Fastest; average_P0_Fast; average_P0_Medium_ISI; average_P0_Slow; average_P0_Slowest];
%Data to be plotted as the error bars
model_error = [disp_P0_Fastest; disp_P0_Fast; disp_P0_Medium_ISI; disp_P0_Slow; disp_P0_Slowest]; 

% Creating axes and the bar graph
if single_or_combined == 1
    figure;
else
    subplot(2,2,4);
end

hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_Exp_cond{i+3}/256;
end
hold off
ax = gca;
ax.YGrid = 'on';

% X and Y labels
ylabel ('Amplitude(µV)');
set(gca,'XTick',[])
% Creating a legend and placing it outside the bar plot
lg = legend('0.25-0.5s','0.5-1s','1-2s','2-4s','4-8s','AutoUpdate','off');
lg.Location = 'BestOutside';
lg.Orientation = 'Horizontal';
hold on;
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none','HandleVisibility','off');
end
errorbar(x(i), model_series(1,1), model_error(1,1), 'r', 'linestyle', 'none','HandleVisibility','off');
title('Amplitudes P0 averaged across dB');
