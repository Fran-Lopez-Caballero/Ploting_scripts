%% Psychoacoustic plots

% I want white backgrounds
set(0,'defaultfigurecolor',[1 1 1]); 

% We will get the data from the Mega_variable, so define which one
gavr_name = 'GAVR_12C_vs_14FE'; % Stats based on subjects from this average
groups_to_plot = {'FE','C'};
channel_data = 'Cz'; % 'Cz', 'cluster' doesn't matter for psychoacoustics
deviation_measure = 'err'; % 'dev' or 'err' for standard dev or error
color_group_sind = [[255 0 0]/256;[0 0 0]/256]; % Specific for SIND
color_group_string = {[255 0 0]/256;[0 0 0]/256}; % Specific for SIND
Quiet_treshold_type = 'Original'; % 'Original' OR 'ChrLab'

% Psychoacoustic vars (may be more in the future)
if strcmp(Quiet_treshold_type,'ChrLab')  
    Quiet_thresholds_L = {'QT_L_125Hz', 'QT_L_250Hz', 'QT_L_500Hz', ...
        'QT_L_1000Hz', 'QT_L_2000Hz','QT_L_4000Hz','QT_L_8000Hz'};
    Quiet_thresholds_R = {'QT_R_125Hz', 'QT_R_250Hz', 'QT_R_500Hz', 'QT_R_1000Hz',...
        'QT_R_2000Hz', 'QT_R_4000Hz','QT_R_8000Hz'};
elseif strcmp(Quiet_treshold_type,'Original')  
    Quiet_thresholds_L = {'QuiT_L_1000', 'QuiT_L_1500', 'QuiT_L_2000', 'QuiT_L_3000', 'QuiT_L_4000'};
    Quiet_thresholds_R = {'QuiT_R_1000', 'QuiT_R_1500', 'QuiT_R_2000', 'QuiT_R_3000', 'QuiT_R_4000'};
end
FD = {'FD_250Hz', 'FD_1000Hz', 'FD_4000Hz'};
ITD = {'ITD_500Hz', 'ITD_1000Hz', 'ITD_2000Hz', 'ITD_4000Hz'};
MD = {'MD_4Hz', 'MD_16Hz', 'MD_64Hz'};
SIND = {'SIND'};
% May add more in the future

% Load mega_variable
load([root_dir '/Statistics/' gavr_name '/Mega_variable_FFR_' channel_data '.mat']);

% QUIET THRESHOLDS BUTTERFLY %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left and right ear, FE and C
figure;
h(1) = subplot(2,2,1); h(2) = subplot(2,2,2); 
h(3) = subplot(2,2,3); h(4) = subplot(2,2,4); 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Mega_variable_FFR(1,:),'Group'));
    pos_group = find(strcmp(Mega_variable_FFR(:,group_col),groups_to_plot{pg}));
    participant = Mega_variable_FFR(pos_group,1);
    % We need a graph
    values_l = [];
    values_r = [];
for p = 1:length(participant)
for qt = 1:length(Quiet_thresholds_L) % Left and right should be the same
        % Position of the variable
        pos_l = find(strcmp(Mega_variable_FFR(1,:),Quiet_thresholds_L{qt}));
        pos_r = find(strcmp(Mega_variable_FFR(1,:),Quiet_thresholds_R{qt}));
        
        % Retrieve values
        if isempty(Mega_variable_FFR{pos_group(p),pos_l})
            values_l(p,qt) = NaN;
        else
            values_l(p,qt) = Mega_variable_FFR{pos_group(p),pos_l};
        end
        if isempty(Mega_variable_FFR{pos_group(p),pos_r})
            values_r(p,qt) = NaN;
        else
            values_r(p,qt) = Mega_variable_FFR{pos_group(p),pos_r};
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
    if strcmp(Quiet_treshold_type,'ChrLab')  
        xticklabels(h(i),{'125','250','500','1000','2000','8000'})
    elseif strcmp(Quiet_treshold_type,'Original')  
        xticklabels(h(i),{'1000','1500','2000','3000','4000'})
    end
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
    group_col = find(strcmp(Mega_variable_FFR(1,:),'Group'));
    pos_group = find(strcmp(Mega_variable_FFR(:,group_col),groups_to_plot{pg}));
    participant = Mega_variable_FFR(pos_group,1);
    % We need a graph
    values_l = [];
    values_r = [];
for p = 1:length(participant)
for qt = 1:length(Quiet_thresholds_L) % Left and right should be the same
        % Position of the variable
        pos_l = find(strcmp(Mega_variable_FFR(1,:),Quiet_thresholds_L{qt}));
        pos_r = find(strcmp(Mega_variable_FFR(1,:),Quiet_thresholds_R{qt}));
        
        % Retrieve values
        if isempty(Mega_variable_FFR{pos_group(p),pos_l})
            values_l(p,qt) = NaN;
        else
            values_l(p,qt) = Mega_variable_FFR{pos_group(p),pos_l};
        end
        if isempty(Mega_variable_FFR{pos_group(p),pos_r})
            values_r(p,qt) = NaN;
        else
            values_r(p,qt) = Mega_variable_FFR{pos_group(p),pos_r};
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
    n(pg) = size(values_l,1); % for legend label
    values_r_average = mean(values_r,1);
    % n(pg) = size(values_r,1); % for legend label
    eval(['n' groups_to_plot{pg} ' = size(values_r,1);'])
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
    if strcmp(Quiet_treshold_type,'ChrLab')  
        errorbar(h(1),[1 2 3 4 5 6 7],values_l_average,dev_l,'color', color_group/256,'HandleVisibility','off')
    elseif strcmp(Quiet_treshold_type,'Original')  
        errorbar(h(1),[1 2 3 4 5],values_l_average,dev_l,'color', color_group/256,'HandleVisibility','off')
    end
    current_title = 'QT Left Ear';
    title(h(1),current_title)
    hold(h(2),'on')
    plot(h(2),values_r_average','-s','MarkerSize',6,'color',color_group/256,'LineWidth', 1.5);
    if strcmp(Quiet_treshold_type,'ChrLab') 
        errorbar(h(2),[1 2 3 4 5 6 7],values_r_average,dev_r,'color', color_group/256,'HandleVisibility','off')
    elseif strcmp(Quiet_treshold_type,'Original') 
        errorbar(h(2),[1 2 3 4 5],values_r_average,dev_r,'color', color_group/256,'HandleVisibility','off')
    end 
    current_title = 'QT Right Ear';
    title(h(2),current_title)
end

% Build legend 
defined_legend = {};
for i = 1:length(participant_group)
    defined_legend{i} = [participant_group{i} ' (n = ' num2str(n(i)) ')'];
end
% Uniform scales, add axes, titles, legend, etc.
for i = 1:2
    hold (h(i),'on')
    ylim(h(i),[floor(current_min),ceil(current_max)]);
    if strcmp(Quiet_treshold_type,'ChrLab') 
        xticklabels(h(i),{'125','250','500','1000','2000','8000'})
    elseif strcmp(Quiet_treshold_type,'Original') 
        xticklabels(h(i),{'1000','1500','2000','3000','4000'})
    end
    if i == 1
        hLeg = legend(h(i),defined_legend);
        ylabel(h(i),'Threshold (dB)');
        xlabel(h(i),'Frequency (Hz)');
    end
end

% FREQ THRESHOLDS BUTTERFLY %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE and C
figure;
h(1) = subplot(1,2,1); h(2) = subplot(1,2,2); 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Mega_variable_FFR(1,:),'Group'));
    pos_group = find(strcmp(Mega_variable_FFR(:,group_col),groups_to_plot{pg}));
    participant = Mega_variable_FFR(pos_group,1);
    % We need a graph
    values = [];
for p = 1:length(participant)
for qt = 1:length(FD)
    % Position of the variable
    pos = find(strcmp(Mega_variable_FFR(1,:),FD{qt}));

    % Retrieve values
    if isempty(Mega_variable_FFR{pos_group(p),pos})
        values(p,qt) = NaN;
    else
        values(p,qt) = Mega_variable_FFR{pos_group(p),pos};
    end        
end   
end
    if pg == 1
        current_max = max(values(:));
        current_min = min(values(:));
    else
        if max(values(:)) > current_max; current_max = max(values(:)); end
        if min(values(:)) < current_min; current_min = min(values(:)); end
    end
    % Plot
    if strcmp(groups_to_plot{pg},'FE')
        hold(h(2),'on')
        plot(h(2),values(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(2),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' FD'];
        title(h(2),current_title)
    elseif strcmp(groups_to_plot{pg},'C')
        hold(h(1),'on')
        plot(h(1),values(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(1),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' FD'];
        title(h(1),current_title)
    end
end
% Uniform scales, add axes, titles, legend, etc.
for i = 1:2
    hold (h(i),'on')
    ylim(h(i),[floor(current_min),ceil(current_max)]);
    NumTicks = 3;
    L = get(h(i),'XLim');
    set(h(i),'XTick',linspace(L(1),L(2),NumTicks))
    xticklabels(h(i),{'250','1000','4000'})
    if i == 1
        ylabel(h(i),'Smallest detectable diff (Hz)')
        xlabel(h(i),'Frequency (Hz)');
    end
end

% FD THRESHOLDS GAVR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE and C
figure; 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Mega_variable_FFR(1,:),'Group'));
    pos_group = find(strcmp(Mega_variable_FFR(:,group_col),groups_to_plot{pg}));
    participant = Mega_variable_FFR(pos_group,1);
    % We need a graph
    values = [];
for p = 1:length(participant)
for qt = 1:length(FD) % Left and right should be the same
        % Position of the variable
        pos = find(strcmp(Mega_variable_FFR(1,:),FD{qt}));
        
        % Retrieve values
        if isempty(Mega_variable_FFR{pos_group(p),pos})
            values(p,qt) = NaN;
        else
            values(p,qt) = Mega_variable_FFR{pos_group(p),pos};
        end
        
end   
end
    if pg == 1
        current_max = max(values(:));
        current_min = min(values(:));
    else
        if max(values(:)) > current_max; current_max = max(values(:)); end
        if min(values(:)) < current_min; current_min = min(values(:)); end
    end
    
    % Correct for NaN
    values = rmmissing(values);
    n(pg) = size(values,1); % For legend
    % Compute GAVR and STDEV/ERR
    values_average = mean(values,1);
    values_stdev = std(values,1);
    values_stderr = values_stdev/(sqrt(size(values,1)));
    if strcmp(deviation_measure,'dev') 
        dev = values_stdev;
    elseif strcmp(deviation_measure,'err')
        dev = values_stderr;
    end

    if strcmp(groups_to_plot{pg},'FE')
        color_group = [255 0 0];
    elseif strcmp(groups_to_plot{pg},'C')
        color_group = [0 0 0];
    end
    plot(values_average','-s','MarkerSize',6,'color',color_group/256,'LineWidth', 1.5);
    hold on
    errorbar([1 2 3],values_average,dev,'color', color_group/256,'HandleVisibility','off')
    current_title = 'Frequency discrimination thresholds';
    title(current_title)
end
% Uniform scales, add axes, titles, legend, etc.
hold on
ylim([floor(current_min),ceil(current_max)]);
NumTicks = 3;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
xticklabels({'250','1000','4000'})
set(gca,'XLim',[0.5 3.5]); 
% Build legend 
defined_legend = {};
for i = 1:length(participant_group)
    defined_legend{i} = [participant_group{i} ' (n = ' num2str(n(i)) ')'];
end
legend(defined_legend);
ylabel('Smallest detectable diff (Hz)')
xlabel('Frequency (Hz)');

% AMP MODULATION THRESHOLDS BUTTERFLY %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE and C
figure;
h(1) = subplot(1,2,1); h(2) = subplot(1,2,2); 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Mega_variable_FFR(1,:),'Group'));
    pos_group = find(strcmp(Mega_variable_FFR(:,group_col),groups_to_plot{pg}));
    participant = Mega_variable_FFR(pos_group,1);
    % We need a graph
    values = [];
for p = 1:length(participant)
for qt = 1:length(MD)
    % Position of the variable
    pos = find(strcmp(Mega_variable_FFR(1,:),MD{qt}));

    % Retrieve values
    if isempty(Mega_variable_FFR{pos_group(p),pos})
        values(p,qt) = NaN;
    else
        values(p,qt) = Mega_variable_FFR{pos_group(p),pos};
    end        
end   
end
    if pg == 1
        current_max = max(values(:));
        current_min = min(values(:));
    else
        if max(values(:)) > current_max; current_max = max(values(:)); end
        if min(values(:)) < current_min; current_min = min(values(:)); end
    end
    % Plot
    if strcmp(groups_to_plot{pg},'FE')
        hold(h(2),'on')
        plot(h(2),values(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(2),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' MD'];
        title(h(2),current_title)
    elseif strcmp(groups_to_plot{pg},'C')
        hold(h(1),'on')
        plot(h(1),values(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(1),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' MD'];
        title(h(1),current_title)
    end
end
% Uniform scales, add axes, titles, legend, etc.
for i = 1:2
    hold (h(i),'on')
    ylim(h(i),[floor(current_min),ceil(current_max)]);
    NumTicks = 3;
    L = get(h(i),'XLim');
    set(h(i),'XTick',linspace(L(1),L(2),NumTicks))
    xticklabels(h(i),{'4','16','64'})
    if i == 1
        ylabel(h(i),'Smallest detectable diff (dB)')
        xlabel(h(i),'Frequency (Hz)');
    end
end

% AMP MODULATION THRESHOLDS GAVR %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE and C
figure; 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Mega_variable_FFR(1,:),'Group'));
    pos_group = find(strcmp(Mega_variable_FFR(:,group_col),groups_to_plot{pg}));
    participant = Mega_variable_FFR(pos_group,1);
    % We need a graph
    values = [];
for p = 1:length(participant)
for qt = 1:length(MD) % Left and right should be the same
        % Position of the variable
        pos = find(strcmp(Mega_variable_FFR(1,:),MD{qt}));
        
        % Retrieve values
        if isempty(Mega_variable_FFR{pos_group(p),pos})
            values(p,qt) = NaN;
        else
            values(p,qt) = Mega_variable_FFR{pos_group(p),pos};
        end
        
end   
end
    if pg == 1
        current_max = max(values(:));
        current_min = min(values(:));
    else
        if max(values(:)) > current_max; current_max = max(values(:)); end
        if min(values(:)) < current_min; current_min = min(values(:)); end
    end
    
    % Correct for NaN
    values = rmmissing(values);
    % Compute GAVR and STDEV/ERR
    n(pg) = size(values,1);
    values_average = mean(values,1);
    values_stdev = std(values,1);
    values_stderr = values_stdev/(sqrt(size(values,1)));
    if strcmp(deviation_measure,'dev') 
        dev = values_stdev;
    elseif strcmp(deviation_measure,'err')
        dev = values_stderr;
    end

    if strcmp(groups_to_plot{pg},'FE')
        color_group = [255 0 0];
    elseif strcmp(groups_to_plot{pg},'C')
        color_group = [0 0 0];
    end
    plot(values_average','-s','MarkerSize',6,'color',color_group/256,'LineWidth', 1.5);
    hold on
    errorbar([1 2 3],values_average,dev,'color', color_group/256,'HandleVisibility','off')
    current_title = 'Amplitude modulation discrimination thresholds';
    title(current_title)
end
% Uniform scales, add axes, titles, legend, etc.
hold on
ylim([floor(current_min),ceil(current_max)]);
NumTicks = 3;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
xticklabels({'4','16','64'})
set(gca,'XLim',[0.5 3.5]); 
% Build legend 
defined_legend = {};
for i = 1:length(participant_group)
    defined_legend{i} = [participant_group{i} ' (n = ' num2str(n(i)) ')'];
end
legend(defined_legend);
ylabel('Smallest detectable diff (dB)')
xlabel('Frequency (Hz)');

% ITD THRESHOLDS BUTTERFLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE and C
figure;
h(1) = subplot(1,2,1); h(2) = subplot(1,2,2); 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Mega_variable_FFR(1,:),'Group'));
    pos_group = find(strcmp(Mega_variable_FFR(:,group_col),groups_to_plot{pg}));
    participant = Mega_variable_FFR(pos_group,1);
    % We need a graph
    values = [];
for p = 1:length(participant)
for qt = 1:length(ITD)
    % Position of the variable
    pos = find(strcmp(Mega_variable_FFR(1,:),ITD{qt}));

    % Retrieve values
    if isempty(Mega_variable_FFR{pos_group(p),pos})
        values(p,qt) = NaN;
    else
        values(p,qt) = Mega_variable_FFR{pos_group(p),pos};
    end        
end   
end
    if pg == 1
        current_max = max(values(:));
        current_min = min(values(:));
    else
        if max(values(:)) > current_max; current_max = max(values(:)); end
        if min(values(:)) < current_min; current_min = min(values(:)); end
    end
    % Plot
    if strcmp(groups_to_plot{pg},'FE')
        hold(h(2),'on')
        plot(h(2),values(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(2),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' ITD'];
        title(h(2),current_title)
    elseif strcmp(groups_to_plot{pg},'C')
        hold(h(1),'on')
        plot(h(1),values(:,:)','-s','MarkerSize',6,'LineWidth', 1.5);
        hLeg = legend(h(1),participant);
        set(hLeg,'visible','off')
        current_title = [groups_to_plot{pg} ' ITD'];
        title(h(1),current_title)
    end
end
% Uniform scales, add axes, titles, legend, etc.
for i = 1:2
    hold (h(i),'on')
    ylim(h(i),[floor(current_min),ceil(current_max)]);
    NumTicks = 4;
    L = get(h(i),'XLim');
    set(h(i),'XTick',linspace(L(1),L(2),NumTicks))
    xticklabels(h(i),{'500','1000','2000','4000'})
    if i == 1
        ylabel(h(i),'Smallest detectable diff (?s)')
        xlabel(h(i),'Frequency (Hz)');
    end
end

% ITD THRESHOLDS GAVR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FE and C
figure; 
% Prepare variable to plot graph
for pg = 1:length(groups_to_plot)
    % Position of subject group
    group_col = find(strcmp(Mega_variable_FFR(1,:),'Group'));
    pos_group = find(strcmp(Mega_variable_FFR(:,group_col),groups_to_plot{pg}));
    participant = Mega_variable_FFR(pos_group,1);
    % We need a graph
    values = [];
for p = 1:length(participant)
for qt = 1:length(ITD) % Left and right should be the same
        % Position of the variable
        pos = find(strcmp(Mega_variable_FFR(1,:),ITD{qt}));
        
        % Retrieve values
        if isempty(Mega_variable_FFR{pos_group(p),pos})
            values(p,qt) = NaN;
        else
            values(p,qt) = Mega_variable_FFR{pos_group(p),pos};
        end
        
end   
end
    if pg == 1
        current_max = max(values(:));
        current_min = min(values(:));
    else
        if max(values(:)) > current_max; current_max = max(values(:)); end
        if min(values(:)) < current_min; current_min = min(values(:)); end
    end
    
    % Correct for NaN
    values = rmmissing(values);
    % Compute GAVR and STDEV/ERR
    n(pg) = size(values,1);
    values_average = mean(values,1);
    values_stdev = std(values,1);
    values_stderr = values_stdev/(sqrt(size(values,1)));
    if strcmp(deviation_measure,'dev') 
        dev = values_stdev;
    elseif strcmp(deviation_measure,'err')
        dev = values_stderr;
    end

    if strcmp(groups_to_plot{pg},'FE')
        color_group = [255 0 0];
    elseif strcmp(groups_to_plot{pg},'C')
        color_group = [0 0 0];
    end
    plot(values_average','-s','MarkerSize',6,'color',color_group/256,'LineWidth', 1.5);
    hold on
    errorbar([1 2 3 4],values_average,dev,'color', color_group/256,'HandleVisibility','off')
    current_title = 'ITD discrimination thresholds';
    title(current_title)
end
% Uniform scales, add axes, titles, legend, etc.
hold on
ylim([floor(current_min),ceil(current_max)]);
NumTicks = 4;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
xticklabels({'500','1000','2000','4000'})
set(gca,'XLim',[0.5 4.5]); 
% Build legend 
defined_legend = {};
for i = 1:length(participant_group)
    defined_legend{i} = [participant_group{i} ' (n = ' num2str(n(i)) ')'];
end
legend(defined_legend);
ylabel('Smallest detectable diff (?s)')
xlabel('Frequency (Hz)');

% SIND BUTTERFLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1 % Silly wat to compress section
pos_measu = find(strcmp(Mega_variable_FFR(1,:),'SIND'));
table_scatter = [];
for pg = 1:length(groups_to_plot)
    group_indices = find(strcmp(Mega_variable_FFR(:,2),groups_to_plot{pg}));
    for i = 1:length(group_indices)
        if isempty(Mega_variable_FFR{group_indices(i),pos_measu})
            table_scatter(i,pg) = NaN;
        else
            table_scatter(i,pg) = Mega_variable_FFR{group_indices(i),pos_measu};
        end
    end
end

% If different numbers of C and FE, it adds zeros to complete tables, correct for that
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(groups_to_plot) ~= 2 % Just to be sure in case we include chronics
    error('If using more than two groups, reprogram next lines');
end
for pg = 1:length(groups_to_plot)
    eval(['group_indices_' num2str(pg) ' = find(strcmp(Mega_variable_FFR(:,2),groups_to_plot{pg}));'])
end
% If they are the same size it won't do anything
if length(group_indices_1) > length(group_indices_2)
    difference =  length(group_indices_1) - length(group_indices_2);
    % So this means first column is 'difference' rows longer than second 
    table_scatter(end+1-difference:end,2) = NaN;
    % Therefore those 'extra' positions at the end of columnn 2 are NaN
elseif length(group_indices_2) > length(group_indices_1)
    difference = length(group_indices_2) - length(group_indices_1);
    % So this means second column is 'difference' rows longer than first 
    table_scatter(end+1-difference:end,1) = NaN;
    % Therefore those 'extra' positions at the end of columnn 1 are NaN
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute an independent sample t test and retrieve p value
[~,p_value,~,t_stats] = ttest2(table_scatter(:,1),table_scatter(:,2));
% Get mean values
mean_plot_left = nanmean(table_scatter(:,1));
mean_plot_right = nanmean(table_scatter(:,2));

% Get std dev (for later use in GAVR)
stddev_plot_left = nanstd(table_scatter(:,1));
stddev_plot_right = nanstd(table_scatter(:,2));

% Get std err (for later use in GAVR)
size1 = sum(~isnan(table_scatter(:,1)),1);
size2 = sum(~isnan(table_scatter(:,2)),1);
stderr_plot_left = stddev_plot_left/(sqrt(size1));
stderr_plot_right = stddev_plot_right/(sqrt(size2));

point_settings = {'MarkerFaceColor',color_group_sind,'MarkerEdgeColor','white','PointSize',80,'LineWidth',1};
plot_settings = [point_settings]; % Something weird about additional wiskers (in case needed)
figure;
[xPositions, yPositions, Label, RangeCut, FigHandles] = UnivarScatter(table_scatter,plot_settings{:});

% set(gcf,'Position',[0,0,600,300])
set(gcf,'Position',[500,250,300,300])
y_title = 'SNR needed for 50% correct words';
y_title = strrep(y_title,'_',' ');
ylabel(y_title); 
xticklabels(groups_to_plot) 
h=gca; h.XAxis.TickLength = [0 0];
h.YGrid = 'on';
h.GridLineStyle = '--';

% Add p value text
if p_value < 0.05
    color_p_value = 'red';
else
    color_p_value = 'black';
end
if p_value < 0.001
    label_pvalue = 'p < 0.001';
else
    str_pvalue = num2str(p_value);
    label_pvalue = ['p = ' str_pvalue(1:5)];
end
title({['\color{' color_p_value '}' label_pvalue '']})

% Add longer mean lines
hold on;
x_values = xlim;
plot([x_values(1)+x_values(1)*0.25 x_values(2)-x_values(2)*0.45],[mean_plot_left,mean_plot_left],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
hold on;
plot([x_values(2)-x_values(2)*0.35 x_values(2)-x_values(2)*0.05],[mean_plot_right,mean_plot_right],'LineWidth',3,'color',[0.5 0.5 0.5 0.5]);
end

% SIND GAVR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 == 1 % silly way to compress into a section
model_series = [mean_plot_left; mean_plot_right];
%Data to be plotted as the error bars
if strcmp(deviation_measure,'dev')
    model_error = [stddev_plot_left; stddev_plot_right]; 
elseif strcmp(deviation_measure,'err')
    model_error = [stderr_plot_left; stderr_plot_right]; 
end

figure;
hold on
for i = 1:length(model_series)
    h=bar(i,model_series(i));
    h.FaceColor = color_group_string{i};
end
set(gca, 'FontSize',12,'XTick',[1 2],'XTickLabel',{['FE (n =  ' num2str(size1) ')'],['C (n = ' num2str(size2) ')']});
ylabel('SNR needed for 50% correct words')
% Finding the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
hold on
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'Color', [0.5 0.5 0.5] , 'linestyle', 'none','HandleVisibility','off');
end
hold off
end