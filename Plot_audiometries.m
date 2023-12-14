%% Psychoacoustic plots

% I want white backgrounds
set(0,'defaultfigurecolor',[1 1 1]); 

groups_to_plot = {'C'};
deviation_measure = 'dev'; % 'dev' or 'err' for standard dev or error

% Psychoacoustic vars (may be more in the future)
Quiet_thresholds_L = {'QT_L_1000Hz', 'QT_L_1500Hz', 'QT_L_2000Hz', ...
    'QT_L_3000Hz', 'QT_L_4000Hz'};
Quiet_thresholds_R = {'QT_R_1000Hz', 'QT_R_1500Hz', 'QT_R_2000Hz', ...
    'QT_R_3000Hz', 'QT_R_4000Hz'};

% Load mega_variable
T1 = readtable('C:/Project/Audiometries/Audiometries_for_stats_Project_v2.xlsx');
table1_columns = T1.Properties.VariableNames;
Audiometry_table = table2cell(T1);
Audiometry_table = [table1_columns; Audiometry_table];

% QUIET THRESHOLDS BUTTERFLY %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left and right ear, FE and C
figure;
h(1) = subplot(2,2,1); h(2) = subplot(2,2,2); 
h(3) = subplot(2,2,3); h(4) = subplot(2,2,4); 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Audiometry_table(1,:),'Group'));
    pos_group = find(strcmp(Audiometry_table(:,group_col),groups_to_plot{pg}));
    participant = Audiometry_table(pos_group,1);
    % We need a graph
    values_l = [];
    values_r = [];
for p = 1:length(participant)
for qt = 1:length(Quiet_thresholds_L) % Left and right should be the same
        % Position of the variable
        pos_l = find(strcmp(Audiometry_table(1,:),Quiet_thresholds_L{qt}));
        pos_r = find(strcmp(Audiometry_table(1,:),Quiet_thresholds_R{qt}));
        
        % Retrieve values
        if isempty(Audiometry_table{pos_group(p),pos_l})
            values_l(p,qt) = NaN;
        else
            values_l(p,qt) = Audiometry_table{pos_group(p),pos_l};
        end
        if isempty(Audiometry_table{pos_group(p),pos_r})
            values_r(p,qt) = NaN;
        else
            values_r(p,qt) = Audiometry_table{pos_group(p),pos_r};
        end
        
end   
end
    if pg == 1
        values = [values_l values_r];
        current_max = max(values(:));
        current_min = min(values(:));
    else
        values = [values_l values_r];
        if max(values(:)) > current_max; current_max = max(values(:)); end
        if min(values(:)) < current_min; current_min = min(values(:)); end
    end
    % See values
    if strcmp(groups_to_plot{pg},'FE')
        hold(h(3),'on')
        plot(h(3),values_l(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(3),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' QT Left Ear'];
        title(h(3),current_title)
        hold(h(4),'on')
        plot(h(4),values_r(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(4),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' QT Right Ear'];
        title(h(4),current_title)
    elseif strcmp(groups_to_plot{pg},'C')
        hold(h(1),'on')
        plot(h(1),values_l(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(1),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' QT Left Ear'];
        title(h(1),current_title)
        hold(h(2),'on')
        plot(h(2),values_r(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(2),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' QT Right Ear'];
        title(h(2),current_title)
    end
%         plot(Grafic_P2(1,:),'color', color_Exp_cond{1}/256, 'LineWidth', 1.5);
%         hold on; plot(values_l(:,:),'-s','MarkerSize',6,'MarkerEdgeColor',[255 0 0]/256,'MarkerFaceColor',[0 0 0]/256,'HandleVisibility','off')
end
% Uniform scales, add axes, titles, legend, etc.
for i = 1:4
    hold (h(i),'on')
    ylim(h(i),[floor(current_min),ceil(current_max)]);
    xticklabels(h(i),{'1000','1500','2000','3000','4000'})
    if i == 3
        ylabel(h(i),'Threshold (dB)')
        xlabel(h(i),'Frequency (Hz)');
    end
end

% QUIET THRESHOLDS GAVR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left and right ear, FE and C
figure;
h(1) = subplot(1,2,1); h(2) = subplot(1,2,2); 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Audiometry_table(1,:),'Group'));
    pos_group = find(strcmp(Audiometry_table(:,group_col),groups_to_plot{pg}));
    participant = Audiometry_table(pos_group,1);
    % We need a graph
    values_l = [];
    values_r = [];
for p = 1:length(participant)
for qt = 1:length(Quiet_thresholds_L) % Left and right should be the same
        % Position of the variable
        pos_l = find(strcmp(Audiometry_table(1,:),Quiet_thresholds_L{qt}));
        pos_r = find(strcmp(Audiometry_table(1,:),Quiet_thresholds_R{qt}));
        
        % Retrieve values
        if isempty(Audiometry_table{pos_group(p),pos_l})
            values_l(p,qt) = NaN;
        else
            values_l(p,qt) = Audiometry_table{pos_group(p),pos_l};
        end
        if isempty(Audiometry_table{pos_group(p),pos_r})
            values_r(p,qt) = NaN;
        else
            values_r(p,qt) = Audiometry_table{pos_group(p),pos_r};
        end
        
end   
end
    if pg == 1
        values = [values_l values_r];
        current_max = max(values(:));
        current_min = min(values(:));
    else
        values = [values_l values_r];
        if max(values(:)) > current_max; current_max = max(values(:)); end
        if min(values(:)) < current_min; current_min = min(values(:)); end
    end
    
    % Correct for NaN
    values_l = rmmissing(values_l);
    values_r = rmmissing(values_r);
    % Compute GAVR and STDEV/ERR
    values_l_average = mean(values_l,1);
    values_r_average = mean(values_r,1);
    values_l_stdev = std(values_l,1);
    values_r_stdev = std(values_r,1);
    values_l_stderr = values_l_stdev/(sqrt(size(values_l,1)));
    values_r_stderr = values_r_stdev/(sqrt(size(values_r,1)));
    if strcmp(deviation_measure,'dev') 
        dev_l = values_l_stdev;
        dev_r = values_r_stdev;
    elseif strcmp(deviation_measure,'err')
        dev_l = values_l_stderr;
        dev_r = values_r_stderr;
    end

    if strcmp(groups_to_plot{pg},'FE')
        color_group = [255 0 0];
    elseif strcmp(groups_to_plot{pg},'C')
        color_group = [0 0 0];
    end
    hold(h(1),'on')
    plot(h(1),values_l_average','-s','MarkerSize',6,'color',color_group/256,'LineWidth', 1.5);
    errorbar(h(1),[1 2 3 4 5],values_l_average,dev_l,'color', color_group/256,'HandleVisibility','off')
    current_title = 'QT Left Ear';
    title(h(1),current_title)
    hold(h(2),'on')
    plot(h(2),values_r_average','-s','MarkerSize',6,'color',color_group/256,'LineWidth', 1.5);
    errorbar(h(2),[1 2 3 4 5],values_r_average,dev_r,'color', color_group/256,'HandleVisibility','off')
    current_title = 'QT Right Ear';
    title(h(2),current_title)
end
% Uniform scales, add axes, titles, legend, etc.
for i = 1:2
    hold (h(i),'on')
    ylim(h(i),[floor(current_min),ceil(current_max)]);
    xticklabels(h(i),{'1000','1500','2000','3000','4000'})
    if i == 1
        hLeg = legend(h(i),groups_to_plot);
        ylabel(h(i),'Threshold (dB)')
        xlabel(h(i),'Frequency (Hz)');
    end
end
